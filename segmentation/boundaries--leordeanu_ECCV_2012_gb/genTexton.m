function map = genTexton(img)

data = load('unitex_6_1_2_1.4_2_64.mat');
lab = rgb2lab(img);
L = lab(:,:,1);
[fim] = fbRun(data.fb,L);
[map] = assignTextons(fim,data.tex);
