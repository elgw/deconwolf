function getmse()

psf = df_readTif('PSF-Bars-stack.tif');

Iref = df_readTif('Bars-stack.tif');
Iref = double(Iref);

I15org = double(df_readTif('Bars-G10-P15-stack.tif'));
I30org = double(df_readTif('Bars-G10-P30-stack.tif'));

I15 = df_readTif('dw_Bars-G10-P15-stack.tif');
I15 = double(I15);
%I15 = I15/0.036920;
I30 = df_readTif('dw_Bars-G10-P30-stack.tif');
I30 = double(I30);
%I30 = I30/0.022979;
L15 = df_readTif('dw_P15_deconvolution_lab.tif');
L15 = double(L15);

L30 = df_readTif('dw_P30_deconvolution_lab.tif');
L30 = double(L30);
whos


[scaling, mseI15org] = fminunc(@(x) mse(Iref, I15org*x), 0.1);
[scaling, mseI30org] = fminunc(@(x) mse(Iref, I30org*x), 0.1);
[scaling, mse15] = fminunc(@(x) mse(Iref, I15*x), 1);
[scaling, mse30] = fminunc(@(x) mse(Iref, I30*x), 1);
[scaling, mseL15] = fminunc(@(x) mse(Iref, L15*x), 1);
[scaling, mseL30] = fminunc(@(x) mse(Iref, L30*x), 1);

%mse15 = mse(Iref, I15);
%mse30 = mse(Iref, I30);
fprintf('mse(org15,ref) = %f\n', mseI15org);
fprintf('mse(org30, ref) = %f\n', mseI30org);
fprintf('mse(P15,ref) = %f\n', mse15);
fprintf('mse(P30, ref) = %f\n', mse30);
fprintf('mse(LP15, ref) = %f\n', mseL15);
fprintf('mse(LP30, ref) = %f\n', mseL30);
%mse30 = mse(Iref, I30);
end


function v = mse(A, B)
 v = sqrt(mean((A(:)-B(:)).^2));
end

