function [gb_thin, gb_fat, texton] = simpleGPb(I)

[gb_thin, ~, gb_fat] = Gb_CSG(I);
texton = genTexton(I);
gb_thin = uint8(round(255*gb_thin));
gb_fat = 255*gb_fat;
texton = uint8(texton);

end