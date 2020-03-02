%GIF conversion
[filename, pathname]=uigetfile('*.avi'); %% selects file for you
mObj = VideoReader([pathname,filename]); 

nFrames = mObj.NumberOfFrames;
vidHeight = mObj.Height;
vidWidth = mObj.Width;

mov(1:nFrames) = struct('cdata', [],'colormap', []);

% Read one frame at a time.
for k = 1 :nFrames
    mov(k).cdata = read(mObj, k);
end

delay1=0.0001;%make this zero if you like it does nothing
name='testing.jpg';%random name, gets overwritten
%cd(pathname);

filename='converted.gif';%you have to change this everytime you run the program

%for loop to make gif
for k = 1:2:nFrames
%to reduce frame rate use    k = 1:2:nFrames etc
    
%this code here is only if you want to resize it    
% imwrite(mov(k).cdata,name);
% I=  imread(name);  
% J = imresize(I, 0.8);
% [f,map] = rgb2ind(J, 16);

%or use this code here to leave it at original size
[f,map] = rgb2ind(mov(k).cdata, 256);%the number here is the number of colors you want

if k == 1
imwrite(f,map,filename,'gif', 'Loopcount',inf,'DelayTime',delay1);
else
imwrite(f,map,filename,'gif','WriteMode','append','DelayTime',delay1);
end
end
