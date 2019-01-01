function R_ans = NormalizeRotationMatrix( R )
    [U,S,V] = svd(R);
    R_ans = U * V';
end