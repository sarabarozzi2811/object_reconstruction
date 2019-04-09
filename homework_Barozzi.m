close all
clear all
clc

%% find edges and lines
%parameter
length_longside=243;

%read the image and change it with some filters 
I = imread('Input_image.jpeg');
figure(1), imshow(I), title('Original image');
I2=localcontrast(I, 0.4, 0.5);

imm = rgb2gray(I2);
Image = imadjust(imm);
y= medfilt2(Image);

%compute 'canny' operator to find edges. The parameter chose come from
%varoius experiments
[edges, th] = edge(imm,'canny');
param = th.*[1.5, 2.5];
    
edges = edge(imm,'canny', param);

figure(2), imshow(edges), title('Edgs obtained with canny operator');

%find the most important lines, parameters are obtained throught different
%experiments
minLineLength_vert = 160;
fillGap = 20;
numPeaks = 200;
NHoodSize = [101 51];
vertical_angles = -90:0.5:89.8;

% find lines 
[H,theta,rho] = hough(edges,'RhoResolution', 1, 'Theta', vertical_angles);
% find peaks in hough transform
P = houghpeaks(H,numPeaks,'threshold',ceil(0.05*max(H(:))), 'NHoodSize', NHoodSize);
% find lines using houghlines
lines_2 = houghlines(edges,theta,rho,P,'FillGap',fillGap,'MinLength',minLineLength_vert);

%display image and lines
figure(3), imshow(y),title('Lines obtained with Hough function'), hold on
max_len = 0;

   for k = 1:length(lines_2)
   xy2 = [lines_2(k).point1; lines_2(k).point2];
   plot(xy2(:,1),xy2(:,2),'LineWidth',2,'Color','green');

   text(xy2(1,1),xy2(1,2), num2str(k), 'Color', 'red')
  
   end
   
figure(11), imshow(y), title('Lines used'), hold on
   for k = [1,2,3,4,5,6,10,19,22,24,29,37,43,47,79]
   xy2 = [lines_2(k).point1; lines_2(k).point2];
   plot(xy2(:,1),xy2(:,2),'LineWidth',2,'Color','green');
   text(xy2(1,1),xy2(1,2), num2str(k),'fontsize',10, 'Color', 'red')
  
   end
   
%% reconstruct shape of horizontal faces

%create vector with initial and ending point of the lines
for i=1: length(lines_2)
 point_a(i,:)= [lines_2(i).point1(1), lines_2(i).point1(2), 1]';
 point_b(i,:)= [lines_2(i).point2(1), lines_2(i).point2(2), 1]';
end

%select line needed to find vanishing point, obtained at the intersection of
%parallel lines

%long parallel line on left horizontal face
l2=cross(point_a(2,:),point_b(2,:));
l1=cross(point_a(1,:),point_b(1,:));

vp_3=cross(l2',l1');
vp_3=vp_3/vp_3(3);

%long parallel line on right horizontal face
l22=cross(point_a(22,:),point_b(22,:));
l4=cross(point_a(4,:),point_b(4,:));

vp_4=cross(l22',l4');
vp_4=vp_4/vp_4(3);

%imaged line at infinite 
l_infinity = cross(vp_3, vp_4)';

% to compute the transformation matrix H, that is an affine trasformation
H_aff = [1 0 0; 0 1 0; l_infinity(1, 1)/l_infinity(1,3) l_infinity(1, 2)/l_infinity(1,3) 1];

%aply transformation to the image
tform = projective2d(H_aff.');
outputImage = imwarp(imm, tform);
figure(4);
imshow(outputImage);
title('Affine trasformation of the given image');

%find perpendicular couples of lines
%on left face
l1=[cross(point_a(1,:),point_b(1,:))]';
l37=[cross(point_a(37,:),point_b(37,:))]';

%on right face
l22=[cross(point_a(22,:),point_b(22,:))]';
l43=[cross(point_a(43,:),point_b(43,:))]';

%vector with all the lines used
linea=[l1,l37,l22,l43];

%normalize and trasform with affine trasformation the lines
for i=1: size(linea,2)
    
    linea_2(:,i)=[linea(1,i)/linea(3,i); linea(2,i)/linea(3,i);linea(3,i)/linea(3,i)];  
    lines(:, i) =  H_aff.' \ linea_2(:,i); 
    %arrey transpose e left division:
    %A.\B is the matrix with elements B(i,j)/A(i,j).
  
end

%extract coefficients from the lines
l11 = lines(1,1);
l12 = lines(2,1);
l21 = lines(1,3);
l22 = lines(2,3);
m11 = lines(1,2);
m12 = lines(2,2);
m21 = lines(1,4);
m22 = lines(2,4);

% Defining variables to find C_starprime
M = [l11*m11 (l11*m12 + l12*m11) ; l21*m21 (l21*m22 + l22*m21)];
b = [-l12*m12;-l22*m22 ];

% X = linsolve(A,B) solves the linear system A*X=B using LU factorization 
%     with partial pivoting when A is square, and QR factorization with 
%     column pivoting otherwise. linsolve warns if A is ill conditioned (for 
%     square matrices) or rank deficient (for rectangular matrices). 
x = linsolve(M,b)

C_starprime=[x(1) x(2) 0,
             x(2)  1   0,
              0    0   0];

% obtaining a suitable rectifying homography
%  svd    Singular value decomposition.
%     [U,S,V] = svd(X) produces a diagonal matrix S, of the same 
%     dimension as X and with nonnegative diagonal elements in
%     decreasing order, and unitary matrices U and V so that
%     X = U*S*V'.
 
[U,S,V] = svd(C_starprime);
H_1 = (U * diag([sqrt(S(1, 1)), sqrt(S(2, 2)), 1]));

%find the euclidian trasformation
H_euc = inv(H_1);

R = rotx(deg2rad(180));

% calculate the composite transformation to rectify the image
% img -> affine -> euclidean -> rotation
H_r = R * H_euc * H_aff;
%apply transformation to the image
tform = projective2d(H_r.');
recstucted_image2 = imwarp(imm, tform);
figure(5);
imshow(recstucted_image2), title('Reconstructed image of horizontal faces');

%% determine relative position of the faces and orientation
%we want to measure the angle between the two faces. 
%find line on left and right face that intersect to have an estimation
%of the angle. 

%long line on the left and right hor face
l1=[cross(point_a(1,:),point_b(1,:))]';
l22=[cross(point_a(22,:),point_b(22,:))]';

%short lines on left and right hor face
l37=[cross(point_a(37,:),point_b(37,:))]';
l43=[cross(point_a(43,:),point_b(43,:))]';

%add trasformation of line 2 that will be used later
l2=[cross(point_a(2,:),point_b(2,:))]';

%vector containig all the line
linea=[l1,l22,l37,l43,l2];
%normalize lines and trasform them with transform them with the composite
%transformation
for i=1: size(linea,2) 
    linea_2(:,i)=[linea(1,i)/linea(3,i); linea(2,i)/linea(3,i);linea(3,i)/linea(3,i)];
    lines_3(:, i) =  H_r.' \ linea_2(:,i); 
end
%estract the first two terms of the first pair of lines
l_1 = lines_3(1:2,1); 
l_2 = lines_3(1:2,2);
% measure the cosine between lines
cos_theta = (l_1.' * l_2)/(norm(l_1,2)*norm(l_2,2));
theta_est = acosd(cos_theta);
    
%estract the first two terms of the second pair of lines
l_3 = lines_3(1:2,3); 
l_4 = lines_3(1:2,4);
% measure the cosine between lines
cos_theta1 = (l_3.' * l_4)/(norm(l_3,2)*norm(l_4,2));
theta_est1 = acosd(cos_theta1);
%estimate final theta doing a media of the two obtained values
theta=(theta_est+theta_est1)/2;

%looking for angle from the border of the image
%line of image border
l_left=[1;0;0];
%select parallel lines on the left face (1 and 2) look for the angle 
%from the border of the image and find the mean value 

%line 1 and line of the border
l_1 = lines_3(1:2,1); 
l_4 = l_left(1:2,:); 
% measure the cosine between lines
cos_theta_left = (l_1.' * l_4)/(norm(l_1,2)*norm(l_4,2));
theta_est_left = acosd(cos_theta_left);
%repeat with the other segment
l_3 = lines_3(1:2,5); 
    
% measure the cosine between lines
cos_theta_left1 = (l_3.' * l_4)/(norm(l_3,2)*norm(l_4,2));
theta_est_left1 = acosd(cos_theta_left1);

theta_from_left=((theta_est_left+theta_est_left1)/2)+180;

%find the vertex of the horizontal right face
%select the corresponding lines and intersect them to obtain the points 
l_r_down=[cross(point_a(47,:),point_b(47,:))]';
l_r_right=[cross(point_a(10,:),point_b(10,:))]';
l_r_up=[cross(point_a(79,:),point_b(79,:))]';
l_r_left=[cross(point_a(24,:),point_b(24,:))]';

lines_face=[l_r_down,l_r_right,l_r_up,l_r_left];

%trasformation of lines with H
for i=1: size(lines_face,2) 
    lines_face(:, i) =  H_r.' \ lines_face(:,i); 
end

x_1 = cross(lines_face(:,2),lines_face(:,3)); %up_right
x_2 = cross(lines_face(:,2),lines_face(:,1)); %down_right
x_3 = cross(lines_face(:,1),lines_face(:,4)); %down_left
x_4 = cross(lines_face(:,4),lines_face(:,3)); %up_left

%normalize
x_1 = x_1 ./ x_1(3,1); 
x_2 = x_2 ./ x_2(3,1);
x_3 = x_3 ./ x_3(3,1);
x_4 = x_4 ./ x_4(3,1);

% length of the longside of horizontal face 
length_longside1 = norm(x_2 - x_1,2);
length_longside2 = norm(x_4 - x_3,2);

% do the average
length_longside_img = (length_longside1 + length_longside2) / 2;

% measure the length of the short side
length_shortside1 = norm(x_4 - x_1,2);
length_shortside2 = norm(x_3 - x_2,2);

% do the average
length_shortside_img = (length_shortside1 + length_shortside2) / 2;

% calculate aspect ratio
aspect_ratio = length_longside_img/length_shortside_img;

%Use two line on the left horizontal face that intersect in the down-left
%vertex of the bandeon to obtain the origin of the left face
l1=[cross(point_a(1,:),point_b(1,:))]';
l29=[cross(point_a(29,:),point_b(29,:))]';

line_L1 =  H_r.' \ l1;
line_L29 =  H_r.' \ l29;
x_origin=cross(line_L1,line_L29);
x_origin = x_origin ./ x_origin(3,1);

%compute the distance between the origin found now and the origin of the
%right horizontal face (x_3)
imaged_distance = norm( x_3-x_origin ,2);

% relative position (in mm) calculated with the known parameter 'length_longside' 
relative_position = imaged_distance * length_longside / length_longside_img;

% relative coordinates (in mm)
relative_coordinates = (x_3-x_origin).* length_longside ./ length_longside_img;

%get the rotation of the right face wrt the image reference frame

% express the relative coordinates in the reference frame of the right face
R_from_img_to_right = rotz(deg2rad(360-theta_from_left)); 
R_last = roty(deg2rad(180));

% relative position of the left face with respect to the right face
% computed as a compination of three rotations
relative_pose_from_right_to_left = R_last*R_from_img_to_right*relative_coordinates;

%% calibration matrix

%use costrain on the homography, on the line at infinity and on the vanishing
%point of the vertical lines.
image_size=max(size(imm));

%find the scaling matrix for the given image
H_scaling = diag([1/image_size, 1/image_size, 1]);

%find l at infinte w.r.t. the scaling factor 
l_inf=l_infinity/l_infinity(1,3);
l_inf = H_scaling.' \ l_inf'; 

%vanishing point of vertical lines
l5=[cross(point_a(5,:),point_b(5,:))];
l6=[cross(point_a(6,:),point_b(6,:))];

vp_v1=cross(l5',l6');
vp_v1=vp_v1/vp_v1(3);

vp_vertical1 = H_scaling * vp_v1;

%matrix to find omega
syms a b c d;
omega = [a 0 b; 0 1 c; b c d];

% scaling factor needed in order to get an homogeneous matrix    
H=H_scaling/H_r;

% first compute the element of l_infinity
l1 = l_inf(1,1);
l2 = l_inf(2,1);
l3 = l_inf(3,1);
   
% vector product matrix
lx = [0 -l3 l2;
      l3 0 -l1;
      -l2 l1 0];
     
h1 = H(:,1);
h2 = H(:,2);
%first and second costrain are the ones on vertical vanishing point and
%line at infinite
% third constraint: h1' w h2 = 0
% fourth equation h1'wh1 = h2' w h2
eqn = [lx(1,:)*omega*vp_vertical1(:,1) == 0;
       lx(2,:)*omega*vp_vertical1(:,1) == 0;
       h1.' * omega * h2 == 0;
       h1.' * omega * h1 == h2.' * omega * h2;];

% equations into matrix form
[A,y] = equationsToMatrix(eqn,[a,b,c,d]);

x=linsolve(A,y);

IAC = double([x(1,1) 0 x(2,1); 0 1 x(3,1); x(2,1) x(3,1) x(4,1)]);

%find parameter for calibration matrix
alfa = sqrt(IAC(1,1));
u0 = -IAC(1,3)/(alfa^2);
v0 = -IAC(2,3);
fy = sqrt(IAC(3,3) - (alfa^2)*(u0^2) - (v0^2));
fx = fy /alfa;
K_de = [fx 0 u0; 0 fy v0; 0 0 1];

% denormalize K
K = H_scaling \ K_de;

% get intrinsic parameter after denormalization
fx = K(1,1);
fy = K(2,2);
u0 = K(1,3);
v0 = K(2,3);
alfa = fx/fy;
 fprintf('The values of the parameter in the K matrix are: fx= %d; fy= %d; u0= %d; v0= %d; alfa= %d;',fx,fy,u0,v0,alfa);
%% localization of camera
%find the rotation of the plane with respect to the camera frame and then
%rotate it wrt the plane. 

%find the coordinate of the vertex of the oriziontal face using the measure
%we know from the text, so we recostruc the shape of the horizontal face
x_dl = [0 0];
x_ul = [0 length_longside];
x_dr = [length_longside/aspect_ratio 0];
x_ur = [length_longside/aspect_ratio length_longside];

%points on the right horizontal face
y_dl=cross(l_r_down,l_r_left);
y_dl=y_dl/y_dl(3);
y_ul=cross(l_r_up,l_r_left);
y_ul=y_ul/y_ul(3);
y_ur=cross(l_r_up,l_r_right);
y_ur=y_ur/y_ur(3);
y_dr=cross(l_r_down,l_r_right);
y_dr=y_dr/y_dr(3);
%we can compute an homography that allows us to obtain the real values
%because we know l-infinity and the aspect
%ration
% %  fitgeotrans Fit geometric transformation to control point pairs.
%       fitgeotrans takes pairs of control points and uses them to infer a
%       geometric transformation.
%    
%       TFORM = fitgeotrans(MOVINGPOINTS,FIXEDPOINTS,TRANSFORMATIONTYPE) fits
%       a geometric transformation to the control point pairs MOVINGPOINTS
%       and FIXEDPOINTS. MOVINGPOINTS is an M-by-2 double matrix containing
%       the X and Y coordinates of control points in the image you want to
%       transform. FIXEDPOINTS is an M-by-2 double matrix containing the X
%       and Y coordinates of control points in the base image.
%       TRANSFORMATIONTYPE can be 'nonreflectivesimilarity', 'similarity',
%       'affine', or 'projective'.
H_omog = fitgeotrans([x_ul; x_dl; x_ur; x_dr], [y_ul(1:2)';y_dl(1:2)'; y_ur(1:2)'; y_dr(1:2)'], 'projective');

H_omog = H_omog.T.'; 
%Homog is the transformation mapping world points to image point.
%we can find the rotation of the plane with respect to the camera with the
%following formulas
h_1=H_omog(:,1);
h_2=H_omog(:,2);
h_3=H_omog(:,3);
lambda = 1 / norm(K \ h_1);
r1 = (K \ h_1) * lambda;
r2 = (K \ h_2) * lambda;
r3 = cross(r1,r2);
% rotation of the world with respect to the camera (R cam -> world)
% the reference in this case is the right face
R = [r1, r2, r3];
%find the nearest orthogonal matrix with the svd operator to correct noises
%in the data
% It is possible to use the SVD of a square matrix A to determine the
% orthogonal matrix O closest to A. The closeness of fit is measured by 
% the Frobenius norm of O ? A. The solution is the product UV?.
% This intuitively makes sense because an orthogonal matrix would have 
% the decomposition UIV? where I is the identity matrix, so that if A = U?V?
% then the product A = UV? amounts to replacing the singular values with
% ones.[from Wikipedia]
[U, ~, V] = svd(R);
R = U * V';
% Compute translation vector. This vector is the position of the plane wrt
% the reference frame of the camera.
T = (K \ (lambda * h_3));

cameraRotation = R.';
%express the position of the camera in the plane rference frame. T is
%expressed in the camera ref frame, so we multiply it for the rotation
cameraPosition = -R.'*T;

%% Display orientation and location of the camera w.r.t. right face

figure(6)
plotCamera('Location', cameraPosition, 'Orientation', cameraRotation.', 'Size', 20);
hold on
pcshow([[x_ul; x_dl; x_ur; x_dr], zeros(size([x_ul; x_dl; x_ur; x_dr],1), 1)], ...
    'red','VerticalAxisDir', 'up', 'MarkerSize', 20);
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Position of the camera w.r.t. right face');

%% localization of camera w.r.t left face
R_from_right_to_left=rotz(deg2rad(theta));
R_camera=R*R_from_right_to_left;

%calculate the traslation of the camera
T_cam_plane_left= T + R*relative_pose_from_right_to_left;

%compute the position of the camera w.r.t the other face
camera_position_left= -R_camera.'*T_cam_plane_left;
camera_orientation=R_camera.';

%plot the camera w.r.t. left, the position is the same as before but we
%change the orientation
figure(7)
plotCamera('Location', cameraPosition, 'Orientation', camera_orientation.', 'Size', 20);
hold on
pcshow([[x_ul; x_dl; x_ur; x_dr], zeros(size([x_ul; x_dl; x_ur; x_dr],1), 1)], ...
    'blue','VerticalAxisDir', 'up', 'MarkerSize', 20);
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Position of the camera w.r.t. left face');
%% Display orientation and location w.r.t. both faces
%compute the shape of the right horizontal face
right_face = [[x_ul; x_dl; x_ur; x_dr], zeros(size([x_ul; x_dl; x_ur; x_dr],1), 1)];

%we want to find the position of the left horizontal face, it is found with a rotation 
left_face = right_face;

left_face = [R_from_right_to_left.'*left_face(1,:).', R_from_right_to_left.'*left_face(2,:).' ...
    , R_from_right_to_left.'*left_face(3,:).', R_from_right_to_left.'*left_face(4,:).'];

%translate the face we found
left_face = left_face + relative_pose_from_right_to_left;
left_face = left_face.';

% plot
figure(8)
plotCamera('Location', cameraPosition, 'Orientation', camera_orientation.', 'Size', 20);
hold on
pcshow(left_face,'red', 'VerticalAxisDir', 'up', 'MarkerSize', 20);
hold on
pcshow(right_face,'blue', 'VerticalAxisDir', 'up', 'MarkerSize', 20);
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Position of the camera w.r.t. both faces');

%% vertical face recostruction
%find the projection matrix and the matrix of the trasformation to obtain
%the recostruction of the right vertical face.
%apply a rotation around y axis to have the image plotted in the right
%position
P = K * [R,T]; 
H_vert_l_sr = inv([P(:,1), P(:,3), P(:,4)]);
R_y = roty(deg2rad(180));
H_vert_l_sr = R_y * H_vert_l_sr;

tform = projective2d(H_vert_l_sr.');
outputImage = imwarp(imm, tform);
figure(9);
imshow(outputImage),
title('Vertical shape reconstruction (right side)');
%find projective matrix for the left vertical face. In this case we use the
%rotation of the camera and the traslation respect to the right to have a
%coherent resutl
P2 = K * [-R_camera,T_cam_plane_left]; % projection matrix
H_vert_r_sr = inv([P2(:,1), P2(:,3), P2(:,4)]);
R_y = roty(deg2rad(180));

% H_vert_sr is the matrix mapping the image to the reconstruction of the
% vertical face
H_vert_r_sr = R_y* H_vert_r_sr;
tform = projective2d(H_vert_r_sr.');
outputImage = imwarp(imm, tform);
figure(10);
imshow(outputImage),
title('Vertical shape reconstruction (left side)');
