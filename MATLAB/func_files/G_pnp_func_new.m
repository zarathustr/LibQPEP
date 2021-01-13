function obj = G_pnp_func_new(in1)
coef_J70 = in1(:,70);
coef_J66 = in1(:,66);
coef_J67 = in1(:,67);
coef_J68 = in1(:,68);
coef_J69 = in1(:,69);
obj = reshape([coef_J66.*-2.0,0.0,-coef_J67,0.0,coef_J68.*-2.0,-coef_J69,-coef_J67,-coef_J69,coef_J70.*-2.0],[3, 3]);
