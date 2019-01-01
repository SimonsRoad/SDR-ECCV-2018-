# SDR-ECCV-2018-
This repository provides the source code of the ECCV 2018 paper --- Semi-Dense Reconstructiong Using a Stereo Event Camera. The sample code aims at helping readers understand the work easily.

## Getting Started

This introduction explains how to configure/run the **matlab + cpp** code.

### Prerequisites

We use the Matrix structures and operations in **Eigen**.

```
http://eigen.tuxfamily.org/index.php?title=Main_Page
```

You need **mexplus** to compile the cpp code.

```
git@github.com:kyamagu/mexplus.git
```

To extract raw data from the rosbag, you need **matlab_rosbag**.

```
https://Joeey@bitbucket.org/Joeey/matlab_rosbag.git
```

In some of the evaluation scripts, **mexopencv** are needed.

```
git@github.com:kyamagu/mexopencv.git
```

For visualization and write images to files, **export_fig** toolbox is needed.

```
https://ch.mathworks.com/matlabcentral/fileexchange/23629-export-fig
```

## Running

To run the code, please follow the steps below.

1. Extract raw data (evetns, raw images, depth map, pose, etc.) from a rosbag. See details in the **matlab_rosbag** repo.

2. Generate stereo observations (rectified event maps and time-surface maps) by running **core**/GenerateObservation_script2_***.m

3. The cpp code is in the **wrapper/**, where EPTAM_mapping.cpp is what you are going to compile and call. The compilation is done in the build.mex. There are three entrances in the EPTAM_mapping.cpp:
```
Entry 1: Single thread version of EstimateInvDepthMap.
Entry 2: Estimate the inverse depth using multiple cores.
Entry 3: Evaluate the objective function at an event coordinate for each hypothesis inverse depth.
```

4. To get 3D reconstruction, run **Test_bash.m**.

5. To get the objective function, run Test_objective2_***.m.

## Authors

* **Yi (Joey) Zhou**
* **cavatinajoey@gmail.com**
