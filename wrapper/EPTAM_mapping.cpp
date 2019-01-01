#include <typeinfo>
#include "mexplus/mxarray.h"
#include <mexplus/dispatch.h>
#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/StdVector>
#include <chrono>
#include <thread>
#include <functional>

using namespace std;
using namespace std::chrono;
using mexplus::MxArray;

#define DEBUG_FLAG false
#define DEBUG_GN false
#define MAX_COST 1e+8
#define MIN_COST 10 // this threshold is to deal with the case in which the residual is close zero because of sparcity of temporal information. 
#define NUM_ITER 10
#define NUM_THREAD 4
#define BRUTE_FORCE_STEP 0.2
#define VAR_EVENT 1e+5
#define SOBEL_SCALE 8
#define VAR_RHO_MAX 0.5

// NOTE: The inverse depth is denoted by d in the code, while is denoted by rho in the paper.

// Back-project 2D pixel to 3D given depth and projection matrix
void cam2World(
	Eigen::Vector2d & x,
	double d,
	Eigen::Matrix<double, 3, 4> & CamP,
	Eigen::Vector3d & p)
{
	double z = 1.0 / d;
	Eigen::Matrix<double,4,1> x_ss;
	x_ss << x(0),
					x(1),
					1,
					1;
	Eigen::Vector4d temp;
	temp << 0, 0, 0, z;
	Eigen::Matrix<double,4,4> P_tilde;
	P_tilde.block<3,4>(0,0) = CamP;
	P_tilde.block<1,4>(3,0) = temp;
  Eigen::Vector4d p_s = z * P_tilde.inverse() * x_ss;
  p = p_s.block<3,1>(0,0) / p_s(3);
}

// Warp a pixel in the RV to one observation (left and right DVS image plane respectively)
void warping2(
	Eigen::Vector2d & x,
	double d,
	Eigen::Matrix<double, 3, 4> & CamP_left,
	Eigen::Matrix<double, 3, 4> & CamP_right,
	Eigen::Matrix<double, 3, 4> & T_left_rv,
	Eigen::Vector2d & x1_s,
	Eigen::Vector2d & x2_s)
{
	// back-project to 3D
	Eigen::Vector3d p_rv;
	cam2World(x, d, CamP_left, p_rv);
	// transfer to left DVS coordinate
	Eigen::Matrix3d R_left_rv = T_left_rv.block<3,3>(0,0);
	Eigen::Vector3d t_left_rv = T_left_rv.block<3,1>(0,3);
	Eigen::Vector3d p_left = R_left_rv * p_rv + t_left_rv;
	// project onto left and right DVS image plane
	Eigen::Vector3d x1_hom = CamP_left.block<3,3>(0,0) * p_left + CamP_left.block<3,1>(0,3);
	Eigen::Vector3d x2_hom = CamP_right.block<3,3>(0,0) * p_left + CamP_right.block<3,1>(0,3);
	x1_s = x1_hom.block<2,1>(0,0) / x1_hom(2);
	x2_s = x2_hom.block<2,1>(0,0) / x2_hom(2);
}

// Check whether the warping results are lying inside the image plane, considering the boundary.
bool InsideImage(
	Eigen::Matrix<double,2,1> & x1_s,
	Eigen::Matrix<double,2,1> & x2_s,
	int width,
	int height,
	int w )
{
	if( x1_s(0) < (w-1)/2 || x1_s(0) > width - (w-1)/2 || x1_s(1) < (w-1)/2 || x1_s(1) > height - (w-1)/2 )
  	return false;
  if( x2_s(0) < (w-1)/2 || x2_s(0) > width - (w-1)/2 || x2_s(1) < (w-1)/2 || x2_s(1) > height - (w-1)/2 )
  	return false;
  return true;
}

// Extract a patch of temporal information from a Exp Decayed SAE. The coordinates of the patch are float, 
// therefore bilinear interpolation is needed. The implementation performs the interpolation on all pixels
// in the patch together as a matrix operation, thus more efficiency than pixel-wise operation.
// param:
// 				img      : the src image.
// 				location : the centre of the patch.
// 				w        : squared patch size.
// 				patch    : return a patch of temporal information.
bool bilinear_patch(
    const Eigen::MatrixXd & img, 
    const Eigen::Vector2d & location, 
    const int w,
    Eigen::MatrixXd & patch )
{
  // compute SrcPatch_UpLeft coordinate and SrcPatch_RightDonw coordinate
	// check patch bourndary is inside img boundary
  Eigen::Vector2i SrcPatch_UpLeft, SrcPatch_RightDonw;
  SrcPatch_UpLeft << floor(location[0]) - (w - 1)/2, floor(location[1]) - (w - 1)/2;
  SrcPatch_RightDonw << floor(location[0]) + (w - 1)/2, floor(location[1]) + (w - 1)/2;
  if(DEBUG_FLAG)
  {
  	cout << "SrcPatch_UpLeft: " << endl << SrcPatch_UpLeft << endl;
  	cout << "SrcPatch_RightDonw: " << endl << SrcPatch_RightDonw << endl;
  }
  if( SrcPatch_UpLeft[0] < 0 || SrcPatch_UpLeft[1] < 0 )
  	return false;
  if( SrcPatch_RightDonw[0] >= img.cols() || SrcPatch_RightDonw[1] >= img.rows() )
  	return false;
  // cout << "--check bourndary (bilinear_patch)" << endl;

  // compute q1 q2 q3 q4
  Eigen::Vector2d double_indices;
  double_indices << location[1], location[0];
  
  std::pair<int,int> lower_indices( floor(double_indices[0]), floor(double_indices[1]) );
  std::pair<int,int> upper_indices( lower_indices.first + 1, lower_indices.second + 1 );
  
  double q1 = upper_indices.second - double_indices[1];
  double q2 = double_indices[1] - lower_indices.second;
  double q3 = upper_indices.first - double_indices[0];
  double q4 = double_indices[0] - lower_indices.first;
  // cout << "--compute q1234" << endl;

  // extract Src patch, size (w+1) * (w+1)
  int w2 = w + 1;
  if(SrcPatch_UpLeft[1] + w >= img.rows() || SrcPatch_UpLeft[0] + w >= img.cols())
  	return false;
  Eigen::MatrixXd SrcPatch = img.block(SrcPatch_UpLeft[1], SrcPatch_UpLeft[0], w2, w2);
  if(DEBUG_FLAG)
  {
  	cout << "SrcPatch: " << endl;
  	cout << SrcPatch << endl;
  }
  // cout << "--extract src patch" << endl;

  // Compute R, size (w+1) * w.
  Eigen::MatrixXd R;
  R = q1 * SrcPatch.block(0,0,w2,w) + q2 * SrcPatch.block(0,1,w2,w);
  // cout << "--Compute R" << endl;

  // Compute F, size w * (w+1).
  patch = q3 * R.block(0,0,w,w) + q4 * R.block(1,0,w,w);
  // cout << "--Compute F" << endl;
  return true;
}

// Import inputs from matlab calls.
void LoadInput(
	const mxArray *prhs[],
	int & numEvents,
	Eigen::MatrixXd & xs,
  int & numObservation,
	int & width, int & height,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vSAEs_left,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vSAEs_right,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vdSAEs_du_left,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vdSAEs_dv_left,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vdSAEs_du_right,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vdSAEs_dv_right,
	vector<Eigen::Matrix<double, 3, 4>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 4> > > & vT_left_rv,
	Eigen::Matrix<double, 3, 4> & CamP_left,
	Eigen::Matrix<double, 3, 4> & CamP_right,
	int & w, string & Lnorm )
{
  // events coordinate xs in the RV
  numEvents = MxArray::to<int>(prhs[0]);
  vector<double> vxs;
	MxArray::to<vector<double> >(prhs[1], &vxs);
	xs = Eigen::Map<Eigen::MatrixXd >(vxs.data(), 2, numEvents);
	
	// number of observation
	numObservation = MxArray::to<int>(prhs[2]);
	
	// width and height
	width = MxArray::to<int>(prhs[3]);
	height = MxArray::to<int>(prhs[4]);

	// a set of decayed temporal map of the left DVS.
	MxArray cell_SAE_left(prhs[5]);
	vSAEs_left.reserve(numObservation);
	for(int i = 0;i < numObservation;i++)
	{
		vector<double> SAE_v = cell_SAE_left.at<vector<double> >(i);
		Eigen::Map<Eigen::MatrixXd> m(SAE_v.data(), height, width);
		vSAEs_left.push_back(m);
	}

	// a set of decayed temporal map of the right DVS.
	MxArray cell_SAE_right(prhs[6]);
	vSAEs_right.reserve(numObservation);
	for(int i = 0;i < numObservation;i++)
	{
		vector<double> SAE_v = cell_SAE_right.at<vector<double> >(i);
		Eigen::Map<Eigen::MatrixXd> m(SAE_v.data(), height, width);
		vSAEs_right.push_back(m);
	}

	// a set of dSAE/du of left DVS
	MxArray cell_dSAE_du_left(prhs[7]);
	vdSAEs_du_left.reserve(numObservation);
	for(int i = 0;i < numObservation;i++)
	{
		vector<double> dSAE_du_v = cell_dSAE_du_left.at<vector<double> >(i);
		Eigen::Map<Eigen::MatrixXd> m(dSAE_du_v.data(), height, width);
		vdSAEs_du_left.push_back(m / SOBEL_SCALE);
	}

	// a set of dSAE/dv of left DVS
	MxArray cell_dSAE_dv_left(prhs[8]);
	vdSAEs_dv_left.reserve(numObservation);
	for(int i = 0;i < numObservation;i++)
	{
		vector<double> dSAE_dv_v = cell_dSAE_dv_left.at<vector<double> >(i);
		Eigen::Map<Eigen::MatrixXd> m(dSAE_dv_v.data(), height, width);
		vdSAEs_dv_left.push_back(m / SOBEL_SCALE);
	}

	// a set of dSAE/du of right DVS
	MxArray cell_dSAE_du_right(prhs[9]);
	vdSAEs_du_right.reserve(numObservation);
	for(int i = 0;i < numObservation;i++)
	{
		vector<double> dSAE_du_v = cell_dSAE_du_right.at<vector<double> >(i);
		Eigen::Map<Eigen::MatrixXd> m(dSAE_du_v.data(), height, width);
		vdSAEs_du_right.push_back(m / SOBEL_SCALE);
	}

	// a set of dSAE/dv of right DVS
	MxArray cell_dSAE_dv_right(prhs[10]);
	vdSAEs_dv_right.reserve(numObservation);
	for(int i = 0;i < numObservation;i++)
	{
		vector<double> dSAE_dv_v = cell_dSAE_dv_right.at<vector<double> >(i);
		Eigen::Map<Eigen::MatrixXd> m(dSAE_dv_v.data(), height, width);
		vdSAEs_dv_right.push_back(m / SOBEL_SCALE);
	}

	// a set of poses of left DVS w.r.t RV
	MxArray cell_T_left_rv(prhs[11]);
	vT_left_rv.reserve(numObservation);
	for(int i = 0;i < numObservation;i++)
	{
		vector<double> T_v = cell_T_left_rv.at<vector<double> >(i);
		Eigen::Map<Eigen::Matrix<double, 3, 4> > T(T_v.data(), 3, 4);
		vT_left_rv.push_back(T);
	}

	// left DVS's projection matrix
	vector<double> vCamP_left;
	MxArray::to<vector<double> >(prhs[12], &vCamP_left);
	 CamP_left = Eigen::Map<Eigen::Matrix<double, 3, 4> >(vCamP_left.data());

	// right DVS's projection matrix
	vector<double> vCamP_right;
	MxArray::to<vector<double> >(prhs[13], &vCamP_right);
	CamP_right = Eigen::Map<Eigen::Matrix<double, 3, 4> >(vCamP_right.data());

	// patch size
	w = MxArray::to<int>(prhs[14]);

	// error metric
	Lnorm = MxArray::to<string>(prhs[15]);
}

// Compute the overall energy of an event coordinate given a hypothesis inverse depth.
double ComputeAggCost(
	Eigen::Matrix<double,2,1> & x,
	double d, int numObservation,
	int width, int height,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vSAEs_left,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vSAEs_right,
	vector<Eigen::Matrix<double, 3, 4>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 4> > > & vT_left_rv,
	Eigen::Matrix<double, 3, 4> & CamP_left,
	Eigen::Matrix<double, 3, 4> & CamP_right,
	int w, string & Lnorm )
{
	double Cr = 0.0;
	int numValid = 0;
	int w2 = w * w;
	//loop through observation
	for(int i = 0; i < numObservation; i++)
	{
		// warping 
		Eigen::Vector2d x1_s, x2_s;
		warping2(x, d, CamP_left, CamP_right, vT_left_rv[i], x1_s, x2_s);

		// check bourndary
		if(!InsideImage(x1_s, x2_s, width, height, w))
		{
			if(DEBUG_FLAG)
				cout << "out of bourndary" << endl;
			continue;
		}

		// extract temporal information
		Eigen::MatrixXd tau1, tau2;
		if( bilinear_patch(vSAEs_left[i], x1_s, w, tau1) && bilinear_patch(vSAEs_right[i], x2_s, w, tau2) )
		{
			//Added by Joey on 1 Mar. This zero check was orginally applied
			//in the Jacobian computation to avoid division over zero. I add it here to generate a nice figure of temporal reisual distribution.
			if( (tau1.array() < 0.001).count() <= 0.75 * w2 && (tau2.array() < 0.001).count() <= 0.75 * w2 )
				numValid++;
		}
		else
		{
			continue;
		}

		if(DEBUG_FLAG)
		{
			cout << "tau1: " << endl;
			cout << tau1 << endl;
			cout << "tau2: " << endl;
			cout << tau2 << endl;
		}

		// compute temporal residual
		if(strcmp(Lnorm.c_str(), "l2") == 0)
			Cr += (tau1 - tau2).squaredNorm();
		else if(strcmp(Lnorm.c_str(), "l1") == 0)
			Cr += (tau1 - tau2).lpNorm<1>();
		else
		{
			cout << "ERROR: Wrong norm is chosen!" << endl;
			exit(-1);
		}
		if(DEBUG_FLAG)
			cout << "Cr: " << Cr << endl;
 	}
 	// cout << "--numValid: " << numValid << endl;
	// if(numValid == 0)
	if(numValid < numObservation / 3)
		Cr = MAX_COST;
	else
		Cr = Cr / numValid;
	return Cr;
}

// Compute \partial tau / \partial d, i.e. the derivatives of the temporal information
// w.r.t the inverse depth d.
void Compute_dtau_dd(
	Eigen::Matrix<double,2,1> & x_r,
	double d,
	Eigen::MatrixXd & dSAE_du,
	Eigen::MatrixXd & dSAE_dv,
	const Eigen::Vector2d & location,
	int w,
	Eigen::Matrix<double, 3, 4> & CamP,
	Eigen::Matrix<double, 3, 4> & T_left_rv,
	Eigen::VectorXd & dtau_dd )
{
	// extract patch from dSAE_du/v
	Eigen::MatrixXd dSAE_du_w, dSAE_dv_w;// w*w
	bool b1 = bilinear_patch(dSAE_du, location, w, dSAE_du_w);
	bool b2 = bilinear_patch(dSAE_dv, location, w, dSAE_dv_w);
	if(DEBUG_FLAG)
		cout << "blinear_patch (Compute_dtau_dd): " << b1 << " " << b2 << endl;
	double w_square = w * w;
	Eigen::Map<Eigen::MatrixXd> temp_col1(dSAE_du_w.data(), w_square, 1);// size w^2 * 1
	Eigen::Map<Eigen::MatrixXd> temp_col2(dSAE_dv_w.data(), w_square, 1);
	Eigen::MatrixXd temp((const int)w_square, 2);// size w^2 * 2
	temp.block(0, 0, w_square, 1) = temp_col1;
	temp.block(0, 1, w_square, 1) = temp_col2;
	if(DEBUG_FLAG)
		cout << "extract patch from dSAE_du/v" << endl;

	// knowns
	double u_r = x_r(0), v_r = x_r(1);
	double p11 = CamP(0,0), p13 = CamP(0,2), p14 = CamP(0,3),
         p22 = CamP(1,1), p23 = CamP(1,2), p24 = CamP(1,3);
  Eigen::Matrix3d R = T_left_rv.block<3,3>(0,0);
  Eigen::Vector3d t = T_left_rv.block<3,1>(0,3);
  double r11 = R(0,0), r12 = R(0,1), r13 = R(0,2),
		     r21 = R(1,0), r22 = R(1,1), r23 = R(1,2),
		     r31 = R(2,0), r32 = R(2,1), r33 = R(2,2);
  double tx = t(0), ty = t(1), tz = t(2);

  // compute du/dd, dv/dd
  double A1 = (p11*r11 + p13*r31)*(u_r - p13)/p11 + (p11*r12 + p13*r32)*(v_r - p23)/p22 + (p11*r13 + p13*r33);
  double B1 = p11*tx + p13*tz + p14;
  double C1 = r31*(u_r - p13)/p11 + r32*(v_r - p23)/p22 + r33;
  double D1 = tz;
  double A2 = (p22*r21 + p23*r31)*(u_r-p13)/p11 + (p22*r22 + p23*r32)*(v_r - p23)/p22 + (p22*r23 + p23*r33);
  double B2 = p22*ty + p23*tz + p24;
  double C2 = C1;
  double D2 = D1;

  double du_dd = (B1*C1 - A1*D1)/pow((C1 + D1 * d),2);
  double dv_dd = (B2*C2 - A2*D2)/pow((C2 + D2 * d),2);
  Eigen::Vector2d temp2;
  temp2 << du_dd,
  				 dv_dd;
  dtau_dd = temp * temp2;
  if(DEBUG_GN)
  {
  	cout << "dSAE_du: " << temp_col1.transpose() << endl;
  	cout << "dSAE_dv: " << temp_col2.transpose() << endl;
  	cout << "du_dd: "   << du_dd << endl;
  	cout << "dv_dd: "   << dv_dd << endl;
  }
}

// compute r = ||tau(x1_s) - tau(x2_s)||_2, and J. Refers to the paper.
bool ComputeResidual_Jacobian(
	Eigen::Matrix<double,2,1> & x,
	double d,
	int width, int height, int w,
	Eigen::MatrixXd & SAE_left,
	Eigen::MatrixXd & SAE_right,
	Eigen::MatrixXd & dSAE_du_left,
	Eigen::MatrixXd & dSAE_dv_left,
	Eigen::MatrixXd & dSAE_du_right,
	Eigen::MatrixXd & dSAE_dv_right,
	Eigen::Matrix<double, 3, 4> & T_left_rv,
	Eigen::Matrix<double, 3, 4> & CamP_left,
	Eigen::Matrix<double, 3, 4> & CamP_right,
	double *r,
	double *J )
{
	// warping 
	double w2 = w * w;
	Eigen::Vector2d x1_s, x2_s;
	warping2(x, d, CamP_left, CamP_right, T_left_rv, x1_s, x2_s);
	// check bourndary
	if(!InsideImage(x1_s, x2_s, width, height, w))
	{
		if(DEBUG_FLAG)
			cout << "out of bourndary" << endl;
		return false;
	}
	// extract temporal information
	Eigen::MatrixXd tau1, tau2;
	if( bilinear_patch(SAE_left, x1_s, w, tau1) && bilinear_patch(SAE_right, x2_s, w, tau2) )// check whether the patch is inside the image
	{
		//check the sparsity of the temporal patch, if too much zero elements, the measurement is not considered.
		if( (tau1.array() < 0.001).count() >= 0.75 * w2 || (tau2.array() < 0.001).count() >= 0.75 * w2 )
			return false;

		// compute the residual
		double residual = sqrt((tau1 - tau2).squaredNorm());
		*r = residual;
		
		// compute jacobian
		Eigen::Map<Eigen::MatrixXd> tau1_vec(tau1.data(), w2, 1);
		Eigen::Map<Eigen::MatrixXd> tau2_vec(tau2.data(), w2, 1);
		Eigen::VectorXd dtau1_dd, dtau2_dd;
		Compute_dtau_dd( x,
										 d,
										 dSAE_du_left,
										 dSAE_dv_left,
										 x1_s,
										 w,
										 CamP_left,
										 T_left_rv,
										 dtau1_dd);
		Compute_dtau_dd( x,
										 d,
										 dSAE_du_right,
										 dSAE_dv_right,
										 x2_s,
										 w,
										 CamP_right,
										 T_left_rv,
										 dtau2_dd);
		Eigen::VectorXd temp1 = tau1_vec - tau2_vec;
		Eigen::VectorXd temp2 = dtau1_dd - dtau2_dd;
		*J = 1.0 / residual * temp1.dot(temp2);
		// debug
		if(DEBUG_GN)
		{
			cout << "residual: " << residual << endl;
			cout << "temp1 * temp2: " << temp1.dot(temp2) << endl;
			cout << "temp2.norm: " << temp2.norm() << endl;
			// cout << "dtau1_dd: " << dtau1_dd.transpose() << endl;
			// cout << "dtau2_dd: " << dtau2_dd.transpose() << endl;
			cout << "J: " << *J << endl;
		}
		return true;
	}
	else
		return false;
}

// Continuously optimize the inverse depth (given d_init) using gauss-newton (with line search).
double NonlinearRefine(
	Eigen::Matrix<double,2,1> & x,//x_r
	double* cost,
	double* var,
	double d_init, double d_boundary,
	int numObservation, int width, int height,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vSAEs_left,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vSAEs_right,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vdSAEs_du_left,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vdSAEs_dv_left,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vdSAEs_du_right,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > & vdSAEs_dv_right,
	vector<Eigen::Matrix<double, 3, 4>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 4> > > & vT_left_rv,
	Eigen::Matrix<double, 3, 4> & CamP_left,
	Eigen::Matrix<double, 3, 4> & CamP_right,
	int w, string & Lnorm )
{
	double d_update = d_init;
	double delta_d = 0.0;

	double d_Lb = d_init - d_boundary;
	double d_Ub = d_init + d_boundary; 

	double Cr_old = ComputeAggCost( x, d_init, numObservation, width, height,
																	vSAEs_left, vSAEs_right, 
																	vT_left_rv,
																	CamP_left, CamP_right,
																	w, Lnorm);
	double Cr_new, nominator, denominator;
	*var = VAR_RHO_MAX;
	*cost = MAX_COST;

	// optimization loop
	for(int i = 0;i < NUM_ITER;i++)
	{
		Cr_new = 0;
		nominator = 0;
		denominator = 0;
		int valid_count = 0;

		// compute J and r at linearization point d_update.
		for(int j = 0;j < numObservation;j++)
		{
			double r = 0;
			double J = 0;
			if( ComputeResidual_Jacobian( 
					x, d_update, width, height, w,
					vSAEs_left[j], vSAEs_right[j],
					vdSAEs_du_left[j],
					vdSAEs_dv_left[j],
					vdSAEs_du_right[j],
					vdSAEs_dv_right[j],
					vT_left_rv[j],
					CamP_left,
					CamP_right,
					&r, &J) )
			{
				nominator += J * r;
				denominator += J * J;
				Cr_new += r * r;
				valid_count++;
			}
			else
				continue;
		}

		// TODO: change this criteria
		if(valid_count < numObservation / 3)
		{
			// no valid observation leads to small information.
			if(DEBUG_GN)
				cout << "terminate 0" << endl;
			return 0.0;//not enough valid observation means the event might be induced by noises.
		}
		else
		{
			Cr_new /= valid_count;
			*var = VAR_EVENT / denominator;//
			*cost = Cr_new;
		}

		// debug
		if(DEBUG_GN)
		{
			cout << "*****************************" << endl;
			cout << "The " << i + 1 << " iteration:" << endl;
			cout << "--Cost: " << Cr_new << endl;
			cout << "--d: " << d_update << endl;
			cout << "--var: " << *var << endl;
		}

		// compute increment
		delta_d = -nominator / denominator;

		// add line search to find optimum step length
		double alpha = 1.0;
		bool bBestStepLength = false;
		for(int t = 0;t < 10;t++)
		{
			double Cr_temp = ComputeAggCost( x, d_update + alpha * delta_d, numObservation, width, height,
																			 vSAEs_left, vSAEs_right, 
																			 vT_left_rv,
																			 CamP_left, CamP_right,
																			 w, Lnorm);
			if(Cr_temp >= Cr_old)
				alpha *= 0.5;
			else
			{
				bBestStepLength = true;
				Cr_new = Cr_temp;
				break;
			}
		}
		
		// not found a good step length, terminate
		if( bBestStepLength == false )
		{
			if(DEBUG_GN)
			 cout << "terminate 3" << endl;
			return d_update;
		}

		// terminate criteria
		if( fabs(Cr_new - Cr_old)/Cr_old < 1e-4 )
		{
			if(DEBUG_GN)
				cout << "terminate 1" << endl;
			return d_update;
		}

		// update d
		double d_temp = d_update + alpha * delta_d;
		if(d_temp > d_Ub || d_temp < d_Lb)//the G-N goes out of boundary, where the optimum is not likely to exist.
		{
			if(DEBUG_GN)
				cout << "terminate 4" << endl;
			return d_update;
		}
		else
		{
			d_update = d_temp;
		}
		Cr_old = Cr_new;

		// debug
		if(DEBUG_GN)
		{
			cout << "--d_update: " << d_update << endl;
			cout << "--delta_d: " << alpha * delta_d << endl;
			cout << "*****************************" << endl;
		}
	}
	// debug
	if(DEBUG_GN)
		cout << "terminate 5" << endl;
	return d_update;
}

// The operator of the Inv Depth estimation, which is called by std::thread for parallel computation.
void EstimateInvDepth(
	int id,
	Eigen::MatrixXd & xs, 
	vector<double> d_s,
	Eigen::VectorXd * d_est,
	Eigen::VectorXd * Cr,
	Eigen::VectorXd * var,
	double d_boundary,
	int numEvents, int numObservation, int width, int height, int w,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vSAEs_left,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vSAEs_right,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vdSAEs_du_left,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vdSAEs_dv_left,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vdSAEs_du_right,
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vdSAEs_dv_right,
	vector<Eigen::Matrix<double, 3, 4>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 4> > > vT_left_rv,
	Eigen::Matrix<double, 3, 4> CamP_left,
	Eigen::Matrix<double, 3, 4> CamP_right,
	string Lnorm )
{
	cout << numEvents << " ----events are distributed to " << id << " thread" << endl;
	int numBruteForceStep = d_s.size();
	auto start = steady_clock::now();
	for(int i = 0; i < numEvents; i++)
  {
  	// the target pixel
  	Eigen::Matrix<double,2,1> x = xs.block<2,1>(0,i);
  	// TODO: This is not very good, think about it!
  	if(x(0,0) < 10 || x(0,0) > width - 10 || x(1,0) < 10 || x(1,1) > height - 10)
  		continue;

  	// 1. find d_init by brute-force searching
  	double d_init = d_s[0];
  	double cost_init = MAX_COST;
  	for(int j = 0; j < numBruteForceStep;j++)
  	{
			// compute cost
  		double cost = ComputeAggCost( x, d_s[j], numObservation, width, height,
																		vSAEs_left, vSAEs_right, 
																		vT_left_rv,
																		CamP_left, CamP_right,
																		w, Lnorm ); 
  		if(cost < cost_init && cost > MIN_COST)
  		{
  			cost_init = cost;
  			d_init = d_s[j];
  		}
  	}
  	if(DEBUG_GN)
  		cout << "found d_init by bruote-force: " << d_init << " ,cost: " << cost_init << endl;; 

  	// This is a bad point, not enough clues for reconstruction.
  	if(cost_init == MAX_COST || cost_init == MIN_COST)
  		continue;

  	// 2. continuous refinement
  	double* variance = new double();
  	double* cr = new double();
  	(*d_est)(i) = NonlinearRefine(x, cr, variance , d_init, d_boundary, numObservation, width, height,
																  vSAEs_left, vSAEs_right,
																  vdSAEs_du_left, vdSAEs_dv_left,
																  vdSAEs_du_right, vdSAEs_dv_right,
																  vT_left_rv,
																  CamP_left, CamP_right,
																  w, Lnorm);
  	(*var)(i) = *variance;
  	(*Cr)(i) = *cr;
  	// (*d_est)(i) = d_init;
  	// cout << "the " << i << " event: " << (*d_est)(i) << endl;
  }
  auto end = steady_clock::now();
  // cout << "----The " << id << " thread uses: " << duration_cast<milliseconds>(end - start).count() / 1000.0 << " seconds\n" ;
}

/******************************************************************************************/
/* Entry 1: Single thread version of EstimateInvDepthMap.
/******************************************************************************************/
MEX_DEFINE(EstimateInvDepthMap)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
  // allocate input
  Eigen::MatrixXd xs;
  int numEvents, numObservation, width, height, w;
  vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vSAEs_left;
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vSAEs_right;
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vdSAEs_du_left;
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vdSAEs_dv_left;
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vdSAEs_du_right;
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vdSAEs_dv_right;
	vector<Eigen::Matrix<double, 3, 4>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 4> > > vT_left_rv;
	Eigen::Matrix<double, 3, 4> CamP_left;
	Eigen::Matrix<double, 3, 4> CamP_right;
	string Lnorm;

  // load input
  auto start = steady_clock::now();
  LoadInput(
  	prhs,
  	numEvents, xs, 
  	numObservation, width, height,
  	vSAEs_left, vSAEs_right,
  	vdSAEs_du_left, vdSAEs_dv_left,
  	vdSAEs_du_right, vdSAEs_dv_right,
  	vT_left_rv,
  	CamP_left, CamP_right,
  	w, Lnorm);
  auto end1 = steady_clock::now();
  if(DEBUG_FLAG)
  {
  	std::cout << "loading input takes: " << duration_cast<milliseconds>(end1-start).count() << " milliseconds\n" ;
  }

  // output some of the inputs
  if(DEBUG_FLAG)
  {
		cout << "numEvents: " << numEvents << endl;
		cout << "x:\n" << xs << endl;
		cout << "numObservation: " << numObservation << endl;
		cout << "width: " << width << endl;
		cout << "height: " << height << endl;
		cout << "vT_left_rv{0}: " << endl;
		cout << vT_left_rv[0] << endl;
		cout << "vCamP_left: " << endl;
		cout << CamP_left << endl;
		cout << "vCamP_right: " << endl;
		cout << CamP_right << endl;
		cout << "patch size: " << w << endl;
		cout << w << endl;
		cout << "Lnorm: " << endl;
		cout << Lnorm << endl;
	}

  //allocate output
	Eigen::MatrixXd DMapList = Eigen::MatrixXd::Zero(5, numEvents);

  // brute-force searching parameters
  double d_Ub = 3;//0.3 m
  double d_Lb = 0.16;//6 m
  double d_boundary = BRUTE_FORCE_STEP;
  vector<double> d_s;
  double d_temp = d_Lb;
  while(d_temp < d_Ub)
  {
  	d_s.push_back(d_temp);
  	d_temp += BRUTE_FORCE_STEP;
  }
  int numBruteForceStep = d_s.size();
  if(DEBUG_FLAG)
  	cout << "numBruteForceStep: " << numBruteForceStep << endl;	

	auto end2 = steady_clock::now();

	// loop
  for(int i = 0; i < numEvents; i++)
  {
  	// the target pixel
  	Eigen::Matrix<double,2,1> x = xs.block<2,1>(0,i);
  	// TODO: This is not good, think about it!
  	if(x(0,0) < 10 || x(0,0) > width - 10 || x(1,0) < 10 || x(1,1) > height - 10)
  		continue;

  	// 1. find d_init by brute-force searching
  	double d_init = d_s[0];
  	double cost_init = MAX_COST;
  	for(int j = 0; j < numBruteForceStep;j++)
  	{
			// compute cost
  		double cost = ComputeAggCost( x, d_s[j], numObservation, width, height,
																		vSAEs_left, vSAEs_right, 
																		vT_left_rv,
																		CamP_left, CamP_right,
																		w, Lnorm ); 
  		if(cost < cost_init && cost > MIN_COST)
  		{
  			cost_init = cost;
  			d_init = d_s[j];
  		}
  	}
  	if(DEBUG_FLAG)
  		cout << "find converging basin" << endl;

  	// this is a bad point, not enough clues for reconstruction
  	if(cost_init == MAX_COST || cost_init == MIN_COST)
  		continue;

  	// 2. continuous refinement
  	double * variance = new double();
  	double * cr = new double();
  	double d_est = NonlinearRefine(x, cr, variance, d_init, d_boundary, numObservation, width, height,
																	 vSAEs_left, vSAEs_right,
																	 vdSAEs_du_left, vdSAEs_dv_left,
																	 vdSAEs_du_right, vdSAEs_dv_right,
																	 vT_left_rv,
																	 CamP_left, CamP_right,
																	 w, Lnorm);
  	cout << "The " << i << " event: " << endl;
  	// double d_est = d_init;
  	
  	if(DEBUG_FLAG)
  		cout << "refined by nonlinear optimization" << endl;
  	// assign the DMap
  	DMapList(0,i) = x(0);
  	DMapList(1,i) = x(1);
  	DMapList(2,i) = d_est;
  	DMapList(3,i) = *variance;
  	DMapList(4,i) = *cr;
  }
  auto end3 = steady_clock::now();
  if(DEBUG_FLAG)
  {
  	cout << "inverse depth estimation takes: " << duration_cast<milliseconds>(end3 - end2).count() << " milliseconds\n" ;
  }

  // reture DMap
  vector<double> result(DMapList.data(), DMapList.data() + DMapList.size());
  plhs[0] = MxArray::from(result);
}

/*********************************************************************************************************/
//Entry 2: Estimate the inverse depth using multiple cores.
/*********************************************************************************************************/
MEX_DEFINE(EstimateInvDepthMap_parallel)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
  // allocate input
  Eigen::MatrixXd xs;
  int numEvents, numObservation, width, height, w;
  vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vSAEs_left;
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vSAEs_right;
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vdSAEs_du_left;
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vdSAEs_dv_left;
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vdSAEs_du_right;
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vdSAEs_dv_right;
	vector<Eigen::Matrix<double, 3, 4>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 4> > > vT_left_rv;
	Eigen::Matrix<double, 3, 4> CamP_left;
	Eigen::Matrix<double, 3, 4> CamP_right;
	string Lnorm;

  // load input
  auto start = steady_clock::now();
  LoadInput(
  	prhs,
  	numEvents, xs, 
  	numObservation, width, height,
  	vSAEs_left, vSAEs_right,
  	vdSAEs_du_left, vdSAEs_dv_left,
  	vdSAEs_du_right, vdSAEs_dv_right,
  	vT_left_rv,
  	CamP_left, CamP_right,
  	w, Lnorm);
  auto end1 = steady_clock::now();
  if(DEBUG_FLAG)
  {
  	std::cout << "loading input takes: " << duration_cast<milliseconds>(end1-start).count() << " milliseconds\n" ;
  }

  // output some of the inputs
  if(DEBUG_FLAG)
  {
		cout << "numEvents: " << numEvents << endl;
		cout << "x:\n" << xs << endl;
		cout << "numObservation: " << numObservation << endl;
		cout << "width: " << width << endl;
		cout << "height: " << height << endl;
		cout << "vT_left_rv{0}: " << endl;
		cout << vT_left_rv[0] << endl;
		cout << "vCamP_left: " << endl;
		cout << CamP_left << endl;
		cout << "vCamP_right: " << endl;
		cout << CamP_right << endl;
		cout << "patch size: " << w << endl;
		cout << w << endl;
		cout << "Lnorm: " << endl;
		cout << Lnorm << endl;
	}

  //allocate output
	Eigen::MatrixXd DMapList = Eigen::MatrixXd::Zero(5, numEvents);

  // brute-force searching parameters
  double d_Ub = 3;//0.3 m
  double d_Lb = 0.16;//6 m
  double d_boundary = BRUTE_FORCE_STEP;
  vector<double> d_s;
  double d_temp = d_Lb;
  while(d_temp < d_Ub)
  {
  	d_s.push_back(d_temp);
  	d_temp += BRUTE_FORCE_STEP;
  }

	auto end2 = steady_clock::now();

	// hyper-thread implementation
	// distribute resources. (Note that the events should be randonly distributed to each thread)
	int numXPerThread = numEvents / NUM_THREAD;
	vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > x_thread;
	vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > d_est_thread;
	vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > cr_thread;
	vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd> > var_thread;
	x_thread.reserve(NUM_THREAD);
	d_est_thread.reserve(NUM_THREAD);
	cr_thread.reserve(NUM_THREAD);
	var_thread.reserve(NUM_THREAD);
	int index = 0;
	if(NUM_THREAD == 1)
	{
		x_thread.push_back(xs);
		d_est_thread.push_back(Eigen::VectorXd::Zero(numEvents));
		cr_thread.push_back(Eigen::VectorXd::Zero(numEvents));
		var_thread.push_back(Eigen::VectorXd::Zero(numEvents));
	}
	else
	{
		for(int i = 0;i < NUM_THREAD;i++)
		{
			x_thread.push_back(Eigen::MatrixXd::Zero(2,numXPerThread));
			for(int j = 0;j < numXPerThread;j++)
			{
				int k = index + rand() % (numEvents - index);
				// swap the index_th column and the k_th column
				Eigen::Matrix<double,2,1> x_temp = xs.block<2,1>(0,index);
				xs.block<2,1>(0,index) = xs.block<2,1>(0,k);
				xs.block<2,1>(0,k) = x_temp;
				x_thread[i].block<2,1>(0,j) = xs.block<2,1>(0,index);
				index++;
			}
			d_est_thread.push_back(Eigen::VectorXd::Zero(numXPerThread));
			cr_thread.push_back(Eigen::VectorXd::Zero(numXPerThread));
			var_thread.push_back(Eigen::VectorXd::Zero(numXPerThread));
		}
	}

	auto end21 = steady_clock::now();
	cout << "Resources distribution takes: " << duration_cast<microseconds>(end21 - end2).count() << " microseconds\n" << endl;

	vector<std::thread> threads;
	for(int i = 0; i < NUM_THREAD;i++)
	{
		threads.emplace_back( std::bind( &EstimateInvDepth, i, x_thread[i], d_s, &d_est_thread[i], &cr_thread[i], &var_thread[i], d_boundary,
																		 numXPerThread, numObservation, width, height, w,
																	   vSAEs_left, vSAEs_right,
																		 vdSAEs_du_left, vdSAEs_dv_left, vdSAEs_du_right, vdSAEs_dv_right,
																		 vT_left_rv,
																		 CamP_left, CamP_right, Lnorm) );
		// cout << "thread " << i << " is created" << endl;
	}

	auto end22 = steady_clock::now();
  cout << "Creating " << NUM_THREAD << " threads takes: " << duration_cast<milliseconds>(end22 - end21).count()/1000.0 << " seconds\n" << endl;

	for( auto& thread : threads )
  {
    if( thread.joinable() )
      thread.join();
  }

  auto end23 = steady_clock::now();
  cout << "Computation takes: " << duration_cast<milliseconds>(end23 - end22).count()/1000.0 << " seconds\n" << endl;

  // draw inverse depth map
  for(int i = 0; i < NUM_THREAD; i++)
  {
  	for(int j = 0; j < numXPerThread; j++)
  	{
  		Eigen::Vector2d x = x_thread[i].block<2,1>(0,j);
  		// cout << "x: " << x << endl;
  		// cout << "d_est_thread: " << (d_est_thread[i])(j);
  		DMapList(0, i * numXPerThread + j) = x(0);
  		DMapList(1, i * numXPerThread + j) = x(1);
  		DMapList(2, i * numXPerThread + j) = (d_est_thread[i])(j);
  		DMapList(3, i * numXPerThread + j) = (var_thread[i])(j);
  		DMapList(4, i * numXPerThread + j) = (cr_thread[i])(j);
  	}
  }

  auto end3 = steady_clock::now();
  if(DEBUG_FLAG)
  {
  	cout << "inverse depth estimation takes: " << duration_cast<milliseconds>(end3 - end2).count() << " milliseconds\n" ;
  }

  // return DMap
  vector<double> result(DMapList.data(), DMapList.data() + DMapList.size());
  plhs[0] = MxArray::from(result);
}

/*********************************************************************************************************/
// Entry 3: Compute the objective function values of an event coordinate at each hypothesis inverse depth.
/*********************************************************************************************************/
MEX_DEFINE(ComputeObjective)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	auto start = steady_clock::now() ;
  // load input
  // event coordinate x in the RV
  double* x_ptr = (double*)mxGetData(prhs[0]);
	Eigen::Matrix<double,2,1> x = Eigen::Map<Eigen::Matrix<double,2,1> >(x_ptr, 2, 1);
	if(DEBUG_FLAG)
		cout << "x:\n" << x << endl;

	// the hypothesis inverse depth(s)
	vector<double> ds = MxArray::to<vector<double> >(prhs[1]);
	if(DEBUG_FLAG)
	{
		cout << "ds: " << endl;
		for(int n = 0; n < ds.size();n++)
			 cout << ds[n] << endl;
	}

	// number of observation
	int numObservation = MxArray::to<int>(prhs[2]);
	if(DEBUG_FLAG)
		cout << "numObservation: " << numObservation << endl;

	// width and height
	int width = MxArray::to<int>(prhs[3]);
	int height = MxArray::to<int>(prhs[4]);
	if(DEBUG_FLAG)
	{
		cout << "width: " << width << endl;
		cout << "height: " << height << endl;
	}
	
	// a set of decayed temporal map of the left DVS.
	MxArray cell_SAE_left(prhs[5]);
	std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vSAEs_left;
	vSAEs_left.reserve(numObservation);
	for(int i = 0;i < numObservation;i++)
	{
		vector<double> SAE_v = cell_SAE_left.at<vector<double> >(i);
		Eigen::Map<Eigen::MatrixXd> m(SAE_v.data(), height, width);
		vSAEs_left.push_back(m);
		// if(DEBUG_FLAG)
		// {
		// 	cout << "vSAEs_left: {" << i << "}" << endl;
		// 	cout << m << endl;
		// }
	}

	// a set of decayed temporal map of the right DVS.
	MxArray cell_SAE_right(prhs[6]);
	std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > vSAEs_right;
	vSAEs_right.reserve(numObservation);
	for(int i = 0;i < numObservation;i++)
	{
		vector<double> SAE_v = cell_SAE_right.at<vector<double> >(i);
		Eigen::Map<Eigen::MatrixXd> m(SAE_v.data(), height, width);
		vSAEs_right.push_back(m);
		// if(DEBUG_FLAG)
		// {
		// 	cout << "vSAEs_right: {" << i << "}" << endl;
		// 	cout << m << endl;
		// }
	}

	// a set of poses of left DVS w.r.t RV
	MxArray cell_T_left_rv(prhs[7]);
	std::vector<Eigen::Matrix<double, 3, 4>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 4> > > vT_left_rv;
	vT_left_rv.reserve(numObservation);
	for(int i = 0;i < numObservation;i++)
	{
		vector<double> T_v = cell_T_left_rv.at<vector<double> >(i);
		Eigen::Map<Eigen::Matrix<double, 3, 4> > T(T_v.data(), 3, 4);
		vT_left_rv.push_back(T);
		if(DEBUG_FLAG)
		{
			cout << "vT_left_rv: {" << i << "}" << endl;
			cout << T << endl;
		}
	}

	// left DVS's projection matrix
	vector<double> vCamP_left;
	MxArray::to<vector<double> >(prhs[8], &vCamP_left);
	Eigen::Matrix<double, 3, 4> CamP_left = Eigen::Map<Eigen::Matrix<double, 3, 4> >(vCamP_left.data());
	if(DEBUG_FLAG)
	{
		cout << "vCamP_left: " << endl;
		cout << CamP_left << endl;
	}

	// right DVS's projection matrix
	vector<double> vCamP_right;
	MxArray::to<vector<double> >(prhs[9], &vCamP_right);
	Eigen::Matrix<double, 3, 4> CamP_right = Eigen::Map<Eigen::Matrix<double, 3, 4> >(vCamP_right.data());
	if(DEBUG_FLAG)
	{
		cout << "vCamP_right: " << endl;
		cout << CamP_right << endl;
	}

	// patch size
	int w = MxArray::to<int>(prhs[10]);
	if(DEBUG_FLAG)
	{
		cout << "patch size: " << w << endl;
		cout << w << endl;
	}

	// error metric
	string Lnorm = MxArray::to<string>(prhs[11]);
	if(DEBUG_FLAG)
	{
		cout << "Lnorm: " << endl;
		cout << Lnorm << endl;
	}
	auto end1 = steady_clock::now() ;

	cout << "data loaded" << endl;

  // Loop: ComputeAggCost
  vector<double> Cr;
  Cr.reserve(ds.size());
  for(int i = 0;i < ds.size();i++)
  {
  	Cr.push_back( ComputeAggCost( x, ds[i], numObservation, width, height,
																  vSAEs_left, vSAEs_right, 
																	vT_left_rv,
																	CamP_left, CamP_right,
																	w, Lnorm ) );
  }
	auto end2 = steady_clock::now() ;
	// if(DEBUG_FLAG)
	// {
		std::cout << "loading input takes: " << duration_cast<milliseconds>(end1-start).count() << " milliseconds\n" ;
		std::cout << "compute cost takes: " << duration_cast<microseconds>(end2-end1).count() << " microseconds\n" ;
	// }

  // reture overall cost
  plhs[0] = MxArray::from(Cr);
}

MEX_DISPATCH