function myshow(image)
figure;
imagesc(abs(image)),colormap gray, axis image,axis off
set(gca,'DataAspectRatio',[1 1 1]); %to make 3d image in org ratio. 
set(gca,'Position',[0 0 1 1]);  %this control the margin size important. !!!
end