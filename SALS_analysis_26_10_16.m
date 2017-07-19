%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Section 1 - Cycle through folders   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start with a folder and get a list of all subfolders.
% Finds and prints names of all PNG, JPG, and TIF images in
% that folder and all of its subfolders.
clc; clear; close all; workspace; format longg; format compact;

%% Plotting options
PlotEllipsoids = 1; % set values = 1 to plot
PlotQuiver = 2;
PlotContour = 1;
PlotEllipse = 1;
PlotStream = 2;
PlotHistogram = 1;
Overlay = 2;
OverlayContour = 1;
mag=2; %set magnification if overlaying on microscopy image
SavePlots = 1; %set = 1 to save all output plots eg. quiver, histograms etc.
         


% Define a starting folder.
start_path = fullfile('F:\Google Drive\MATLAB\Stress Relaxation\SRBC2');
% Ask user to confirm or change.
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
    return;
end

% Get list of all subfolders.
allSubFolders = genpath(topLevelFolder);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};

% Loop true folders
while true
    [singleSubFolder, remain] = strtok(remain, ';');
    if isempty(singleSubFolder)
        break;
    end
    listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames)

%% Process all image files in those folders.
for j = 1 :numberOfFolders
    
    %%change current directory (1 at a time) to run code in
    cd(listOfFolderNames{j})
    pwd
    
    if exist('SALS.1.png', 'file')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%   Section 2 - Setup dimensions/resolution  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Method 1 - Popup Dialog box to enter sample details %
        % prompt = 'Number of images? ';
        % N = input(prompt)
        % prompt = 'Sample width (mm)? ';
        % sample_width = input(prompt)
        % prompt = 'Sample height (mm)? ';
        % sample_height = input(prompt)
        % prompt = 'resolution (um)? ';
        % resolution = input(prompt)
        % resolution = resolution/1000
        
        % Method 2 - loop with commonly used numbers/dimensions/resolutions %
        % Temporary method to cycle through the right number of files and determine
        % sample dimensions based on commonly used interogation regions/resolutions
        % give a range for num_images to account for graph.png files etc.
        files = dir('*.png');
        num_images = length(files);
        
                if num_images > 259 && num_images < 400
                    N=260;
                    sample_width= 1.5;
                    sample_height= 0.975;
                    resolution = 0.075
                elseif num_images > 191 && num_images < 210
                    N=192;
                    sample_width = 2;
                    sample_height = 1.5;
                    resolution = 0.125
                elseif num_images > 95 && num_images < 190
                    N=96;
                    sample_width = 1.5;
                    sample_height = 1;
                    resolution = 0.125
                elseif num_images > 63 && num_images < 74
                    N=64;
                    sample_width = 2;
                    sample_height = 2;
                    resolution = 0.25
                elseif num_images > 23 && num_images < 48
                    N=24;
                    sample_width = 1.5;
                    sample_height = 1;
                    resolution = 0.25
                elseif num_images > 15 && num_images < 21
                    N=16;
                    sample_width = 1;
                    sample_height = 1;
                    resolution = 0.25
                else
                    fprintf('\nNumber of images does not relate to previously defined sample width/height/resolution\n')
                end
        
        % Method 3 - Hard Code sample numbers/dimensions/resolutions
%         N=3500;
%         sample_width=2.1;
%         sample_height=1.5;
%         resolution = 0.030;
        
        %% Initialise values/Create arrays/Set file name to be incremented. Increases performance/speed%%
        IMAGE = cell(1,N);
        
        maxints = zeros(1,N);
        angle = zeros(1,N);
        minints =zeros(1,N);
        minangle = zeros(1,N);
        ratio = zeros(1,N);
        maxints2 = zeros(1,N);
        minints2 = zeros(1,N);
        profile1 = cell(360,1);
        circumferenceX = zeros(360,1);
        circumferenceY = zeros(360,1);
        coordinatesX = zeros(360,1);
        coordinatesY = zeros(360,1);
        ints = zeros(1,360);
        Fibreangle = zeros(1,N);
        theta = zeros(1,360);
        xx = zeros(1,N);
        yy = zeros(1,N);
        xx2 = zeros(1,N);
        yy2 = zeros(1,N);
        u1 = zeros(1,N);
        v1 = zeros(1,N);
        ii=1;
        Xstore = cell(1,N);
        Ystore = cell(1,N);
        Zstore = cell(1,N);
        stats(N) = struct('Area', 0, 'MajorAxisLength', 0,'MinorAxisLength', 0, 'Eccentricity', 0, 'Orientation', 0);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Section 3 - Read image & image processing  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Image prefix to read in %
        FileFMT = 'SALS.%d.png';
       
        % N=1 %Debugging single photo       
        clear wait
        wait = waitbar(0,'Please wait...');
        % Cycle through images using sprintf command (1:N eg. file.1, file.2 etc.)
%         N=8
        for i=1:N
            waitbar(i/N)
            IMAGE{i} = imread(sprintf(FileFMT,i)); %Read in image
            I=IMAGE{i};
            
            Icrop = imcrop(I,[200 0 800 720]); %Crop image to remove refections from enclosure
            Igreen=Icrop(:,:,2);    %Consider one channel (converts to greyscale). Testing showed green channel best
            I2 = fliplr(Igreen); %Image taken from reverse side so flip to get angle when looking from laser side
            
%             figure;          %Debugging
%             imshow(I2);
            
            GausI2 = imgaussfilt(I2,4); %apply a gaussian filter to remove noise and smooth the greyscale image
            
%                         figure;          %Debugging
%                         imshow(GausI2);
            
            % Thresholding/processing functions
            greythresh = multithresh(GausI2,2); %Automated multithresh to seperate black, red and 'white' light regions of image
            A=double(greythresh(2)); %Use the second thresh to only account for brightest light region
            binarythresh=(1/256)*A; %Thresh needs to be converted to 0-1 range
            
            B1 = im2bw(GausI2,binarythresh); %binary image
            
%                         figure;          %Debugging
%                         imshow(B1);
            
            %	se = strel('disk', 4); %Set operator shape and size for Open/Close function (smooths image)
            %	B1 = imopen(B1,se);
            
%             figure;          % Debugging
%             imshow(B1);
            
            cleanB1 = bwareaopen(B1, 2000); % Remove blobs smaller than set size from image
            
            stats(i) = regionprops(cleanB1, 'Area','Eccentricity','Orientation','MajorAxisLength', 'MinorAxisLength'); % Blob stats - area used for removing 'non fibres' in plot at end. The rest is just for reference
%             h= figure;          % Debugging
%             imshow(cleanB1);
            
            %% Blob analysis %%
            [labeledImage, numberOfBlobs] = bwlabel(cleanB1, 8); % Label blobs for processing
            blobMeasurements = regionprops(labeledImage, 'BoundingBox'); % Get smallest bounding box for the blob
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%   Section 4 - Determine interogation region & find intensity/angle  %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Display image for bounding box at later step with 'hold on' (Uncomment 'plot' function too to display)
            % figure;           % Debugging
            % imshow(cleanB1, []);
            % title('Selected Image with Bounding Boxes');
            % hold on;
            
            % Loop through all blobs, creating Bounding Box (k should = 1 ie. 1 blob(scatter))
            for k = 1 : numberOfBlobs
                boundingBox = blobMeasurements(k).BoundingBox;	% Get box.
                x1 = boundingBox(1); % boundingBox(1/2) find first corner
                y1 = boundingBox(2);
                x2 = x1 + boundingBox(3) - 1; % boundingBox(3/4) find width/height
                y2 = y1 + boundingBox(4) - 1;
                verticesX = [x1 x2 x2 x1 x1];
                verticesY = [y1 y1 y2 y2 y1];
                
                %plot(verticesX, verticesY, 'g-', 'LineWidth', 2); % Plot the box in the overlay (Debugging)
            end
            
            %% Preparing/Finding intensity distribution about 360 degrees
            % Circle radius for cycling through scattered data determined by plot/bounding box size
            radius = sqrt((boundingBox(3)/2)^2 + (boundingBox(4)/2)^2); % Max distance from centre to box perimeter
            
            [x,y] = find(cleanB1); % Get all non zero pixel x/y locations
            CenterOfMassXY = [mean(x) mean(y)] ; % Centre point found from average of pixel locations
            CenX=CenterOfMassXY(1,2);
            CenY=CenterOfMassXY(1,1);
            
            % Determine distance from centre of circle to outer circumference coordinates of each segment
            for n=1:360
                circumferenceX(n,1)= radius*cos(n*(2*3.14)/360);
                circumferenceY(n,1)= (-1)*radius*sin(n*(2*3.14)/360);
            end
            
            % Global coordinates of circumfence points found by adding on distance from global origin to centre of circle
            for n=1:360
                coordinatesX(n,1)= circumferenceX(n,1)+(CenX);
                coordinatesY(n,1)= circumferenceY(n,1)+(CenY);
            end
            
            % Create vectors Xi,Yi of equal lengths by specifying line end points (circle circumference points)
            for n=1:360
                Xi=[CenX,coordinatesX(n,1)];
                Yi=[CenY,coordinatesY(n,1)];
                
                % Plot the intensity profile of each line created using Xi & Yi (B1 or I1 for binary or greyscale)
                profile1{n,1} = improfile(cleanB1,Xi,Yi);
                
                %         hold on;
                %         plot([570,403],[395,403],'Color','g','LineWidth',1)
                % Determine the average of each profile line created from 0-360 degrees
                ints(n)=mean(profile1{n,1});
            end
            
            %% Plot previous average intensity from 0-360 degrees with a Y axis from 0-255(greyscale)(Change cleanB1 to I2 above) or 0-1(binary)
%                  figure;
%                  plot(ints);
%                  axis([0 360 0 1])
%                  title('Intensity vs Angle data');
%                  xlabel('Angle, ?(°)'); ylabel('Intensity');
            
            % Find max intensity and corresponding angle from 0-360 degrees in an image
            [maxints(ii), angle(ii)] = max(ints);
            [minints(ii), minangle(ii)] = min(ints);
            ratio(ii) = minints(ii)/maxints(ii);
            
            maxints2(ii) = maxints(ii); % generate a non-zero copy to be used in contour plot
            minints2(ii) = minints(ii);
            
            %%%%% CORRECTION FACTOR %%%%%%%%%%%%
            %if sample is not perfectly orientated on slide, create file called in same
            %directory 'correction.txt' specifcy what angle to adjust for. Code will
            %search automatically if file exists within folder
            if exist([pwd filesep 'correction.txt'], 'file')==2
                
                fileID = fopen('correction.txt','r');
                correction = fscanf(fileID, '%d');
                fclose(fileID);
                
                angle(ii)=angle(ii)+correction;
            end
            ii=ii+1; % Increment ii for subsequent image results to be stored
            
            % Scatterplateangle = angle of highest intensity peak in each image
            scatterplateangle=(angle);
            
            % dominantscatterplateangle = average angle over all images (Avg angle of tested sample)
            dominantscatterplateangle=mean(angle);
            
            %% Transform scatter to ??180 degree region
            for l=1:N
                if angle(l)>180
                    angle(l)=angle(l)-180;
                end
            end
            
            % subtract 90 degrees to determine fibre angle between +90 and -90 deg
            for m=1:N
                Fibreangle(m)=angle(m)-90;
            end
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%      Section 5 - Summarising      %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            FibreanglePos = Fibreangle(Fibreangle > 0);
            FibreangleNeg = Fibreangle(Fibreangle < 0);
            FibreangleZero = Fibreangle(Fibreangle == 0);
            
            NumberPos = length(FibreanglePos);
            NumberNeg = length(FibreangleNeg);
            NumberZero = length(FibreangleZero);
            
            % Intensity = max intensity of each image
            Intensity=(maxints);
            % Eccentricity based on Haskett paper. 1 = straight line, 0 = circle
            Eccentricity2 = (2*sqrt(((maxints./2).^2)-((minints./2).^2)))./maxints; 
            
        end
        close(wait)
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%        Section 6 - Plotting       %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         resolution=0.25; %RESET dimensions for testing/Debugging
        %         sample_width=2.25;
        %         sample_height=1.5;
        
        % find u and v for vectors based on angles and eccentricity
        Radian_Fibreangle = Fibreangle*(pi/180);
        for i=1:N
            theta(i) = Radian_Fibreangle(i);
            %r=maxints(i)/minints(i); % scale unit vector components by alignment of SALS scatter
            r=Eccentricity2(i);
            u1(i) = r*cos(theta(i)); % convert polar (theta,r) to cartesian
            v1(i) = r*sin(theta(i));
            
            % choosing to remove fibres below a threshold based on
            % eccentricity and area
            if Eccentricity2(i) < 0.3 && stats(i).Area < 19000
                %Could also base 'fibre elimination' on intensity of input
                %images
                u1(i) = 0;
                v1(i) = 0;
                
                maxints(i) = 0;
                minints(i) = 0;
            end
        end
        
        %Set up
        col = sample_width/resolution;
        row = sample_height/resolution;
        
        %Origin basis for ellipsoid and quiver plots
        A = [resolution:resolution:sample_width];
        A = repmat(A, 1, row);
        B = [resolution:resolution:sample_height];
        B = repmat(B, col, 1);
        B = reshape(B, N, 1);
        
        %%
        %%%%  ELLIPSOID  %%%%
        tensor=zeros(N,6);
        tensor(:,1)=A;
        tensor(:,2)=B;
        tensor(:,4)=maxints;
        tensor(:,5)=minints;
        
        %Summarising information for ellipsoid plot
        A2(:,1)=abs(tensor(:,1));
        B2(:,1)=abs(tensor(:,2));
        C2(:,1)=abs(tensor(:,3));
        u2(:,1)=(abs(tensor(:,4)))*resolution*0.5; %scale size of ellipsoids according to resolution and sample width
        v2(:,1)=(abs(tensor(:,5)))*resolution*0.5;
        w2(:,1)=0.01+abs(tensor(:,6));
        
        
        if PlotEllipsoids == 1
            for s=1:N
                [X,Y,Z]=ellipsoid(0,0,0,u2(s),v2(s),w2(s));
                
                h=surf(X,Y,Z);
                %figure
                rotate(h,[0 0 1],-Fibreangle(s)) %rotate single ellipses by fibre angle, done individually to not affect to every other ellipsoid. negative fibre angle as y axis is reverese to match laser scanning route
                
                XData=get(h,'XData'); %Get X/Y/Z data for individual ellipsoids to make main plot later
                YData=get(h,'YData');
                ZData=get(h,'ZData');
                Xstore{s}=XData+A2(s); %Store above + add on X/Y/Z origin for main plot
                Ystore{s}=YData+B2(s);
                Zstore{s}=ZData+C2(s);
                
                close all; %close individual plots
            end
            for s=1:N
                hold on
                h(s)=surf(Xstore{s},Ystore{s},Zstore{s}); % main plot with all ellipsoids from above
                hold off
            end
            
            view([0 90]);
            set(gca,'GridLineStyle','none')
            shading interp
            colormap([1 0 0])
            lighting phong
            light('Position',[0 0 1],'Style','infinite','Color',[1 1 1]);
            xlabel('Sample width (mm)') % x-axis label
            ylabel('Sample height (mm)') % y-axis label
            axis([0,sample_width+resolution,0,sample_height+resolution])
            set(gca,'Ydir','reverse')
            if SavePlots == 1
                saveas(gcf,'ellipsoid.fig')
                saveas(gcf,'ellipsoid.png')
            else
                disp('save images not selected');
            end
        else
            disp('Ellipsoid plot not selected');
        end
        
        %%
        %%%%  Histology + Quiver  %%%%
        if Overlay == 1
            if mag == 2
                overlay_scaling = 0.59*1000;  %2x Multiply by 1000 to change from mm to microns
            elseif mag == 4
                overlay_scaling = 1.18*1000;  %4x
            elseif mag == 10
                overlay_scaling = 5.9*1000;   %10x
            elseif mag == 20
                overlay_scaling = 11.8*1000;  %20x
            elseif mag == 40
                overlay_scaling = 23.6*1000;  %40x
            else
                disp('set mag/scaling factor');
            end
            
            resolution=resolution*1000; %convert to microns
            sample_width=sample_width*1000;
            sample_height=sample_height*1000;
        else
            overlay_scaling = 1;
        end
        
        A3=A*overlay_scaling; %micron to pixel scaling factor for 2x mag
        B3=B'; % convert previous B to single column vector
        B3=B3*overlay_scaling; % automatically scaled if overlaying
        
        A4=A2'*1000; % mm to microns used for contour overlay. A4 used so not to change previous A values
        B4=B2'*1000;
        
        if PlotQuiver == 1
            if Overlay == 1
                I1=imread('PLM1.tif');
                I2=imread('PLM2.tif');
                I=imadd(I1,I2);
                I=rot90(I,2);
                figure(2);imshow(I);
                
                hold on;
                
                I=imread('Full.tif');
                I=rot90(I,2); %flips image left-right, up-down to correct for microscope
                figure(3);h1=imshow(I);
                
                hold on;
            else
                disp('no overlay for quiver selected');
            end
            %             test = 0.001; % correcting for scaling issues
            %             test2 = 1.09;
            test = 1;
            test2 = 1;
            figure(2)
            h1 = quiver(A3*test2,B3*test2,u1*overlay_scaling*resolution*0.9*test,-v1*overlay_scaling*resolution*0.9*test,0,'.'); %v1 is negative as y axis is reversed below to match SALS scanning route
            c = h1.Color;
            if Overlay == 1
                h1.Color = 'white';
            else
                h1.Color = 'red';
                xlabel('Sample width') % x-axis label
                ylabel('Sample height') % y-axis label
            end
            set(h1,'linewidth',2);
            axis([0,sample_width+resolution,0,sample_height+resolution])
            if Overlay == 1
                axis tight
            end
            set(gca,'Ydir','reverse')
            axis equal
%             if SavePlots == 1
%                 saveas(gcf,'Quiver_Overlay1.fig')
%                 saveas(gcf,'Quiver_Overlay1.png')
%             else
%                 disp('save images not selected');
%             end
            hold off;
            
            if Overlay == 1
                figure(3)
                h1 = quiver(A3*test2,B3*test2,u1*overlay_scaling*resolution*0.9*test,-v1*overlay_scaling*resolution*0.9*test,0,'.');%v1 is negative as y axis is reversed below to match SALS scanning route
                c = h1.Color;
                h1.Color = 'white';
            end
            set(h1,'linewidth',2);
            axis([0,sample_width+resolution,0,sample_height+resolution])
            if Overlay == 1
                axis tight
            end
            
            hold off;
            set(gca,'Ydir','reverse')
            axis equal
%             if SavePlots == 1
%                 saveas(gcf,'Quiver_Overlay2.fig')
%                 saveas(gcf,'Quiver_Overlay2.png')
%             else
%                 disp('save images not selected');
%             end
            
        else
            disp('Quiver plot not selected');
        end
        
        %%
        %%%%  Contour + Quiver  %%%%
        if PlotContour == 1
            %Create x,y,z values for colour contour plot that can be overlaid on quiver
            [xcontour, ycontour] = meshgrid(0:resolution:sample_width+resolution,0:resolution:sample_height+resolution);
            
            % %Select either angle or eccentricity for contour plot Z
            Eccentricity1 = minints2./maxints2;
            Eccentricity2 = (2*sqrt(((maxints2./2).^2)-((minints2./2).^2)))./maxints2;
            
            %Use Eccentricity2 or Fibreangle below for Zcontour
            Zcontour = myvec2mat(Eccentricity2, row, col); %Colour based on Fibreangle. Needs to be in matrix like x and y contour
            paddedZcontour = padarray(Zcontour,[1 1],'replicate'); %add single row/col boarder to matrix so quiver fits inside contour plot
            
            figure(4);
            [cv, ch] = contourf(xcontour, ycontour, paddedZcontour, 300); %number of levels
            set(ch,'edgecolor','none'); %remove gradient/contour lines for smoother image
            hold on; %for overlaying quiver
            
            if OverlayContour == 1
                figure(4)
                h1 = quiver(A3,B3,u1*resolution*0.9,-v1*resolution*0.9,0,'.'); %v1 is negative as y axis is reversed below to match SALS scanning route.Might need to change to A4/B4
                c = h1.Color;
                h1.Color = 'white';
                xlabel('Sample width') % x-axis label
                ylabel('Sample height') % y-axis label
                set(h1,'linewidth',1.5);
                axis([0,sample_width+resolution,0,sample_height+resolution]);
                axis equal
            end
            
            hold off;
            set(gca,'Ydir','reverse')
            c = colorbar;
            c.Label.String = 'Eccentricity';
            caxis([0.45 0.9])
            axis equal
            if SavePlots == 1
                saveas(gcf,'Contour.fig')
                saveas(gcf,'Contour.png')
            else
                disp('save images not selected');
            end
            
        else
            disp('No quiver for contour selected');
        end
        
        %%
        %%%  OPEN ELLIPSE %%%%
        if PlotEllipse == 1
            if Overlay == 1
                I1=imread('PLM1.tif');
                I2=imread('PLM2.tif');
                I=imadd(I1,I2);
                I=flipud(I);
                I=rot90(I,2);
                figure(5);imshow(I);
                hold on;
                
                centerX=A*overlay_scaling;
                centerY=B*overlay_scaling;
                
                majorAxis = maxints*resolution*0.5; %scale size of ellipse based on resolution so no overlap
                minorAxis = minints*resolution*0.5;
            else
                centerX=A*1000;
                centerY=B*1000;
                resolution=resolution*1000;
                figure(5);
                
                majorAxis = maxints*resolution*0.9; %scale size of ellipse based on resolution so no overlap
                minorAxis = minints*resolution*0.9;
            end
            
            theta = linspace(0, 2*pi, 150);
            orientation=-(Fibreangle*pi/180); %angle is negative as y axis is reversed below to match SALS scanning route
            xx=zeros(N,150);
            
            for i=1:N
                xx(i,:) = (majorAxis(i)/2) * sin(theta) + centerX(i);
                yy(i,:) = (minorAxis(i)/2) * cos(theta) + centerY(i);
            end
            
            for i=1:N
                xx2(i,:) = (xx(i,:)-centerX(i))*cos(orientation(i)) - (yy(i,:)-centerY(i))*sin(orientation(i)) + centerX(i);
                yy2(i,:) = (xx(i,:)-centerX(i))*sin(orientation(i)) + (yy(i,:)-centerY(i))*cos(orientation(i)) + centerY(i);
                
                h3=plot(xx2(i,:),yy2(i,:));
                c1 = h3.Color;
                if Overlay == 1
                    h3.Color = 'white';
                else
                    h3.Color = 'red';
                end
                set(h3,'linewidth',1.25);
                hold on;
            end
            set(gca,'Ydir','reverse')
            axis equal
            if SavePlots == 1
                saveas(gcf,'Open_Ellipse.fig')
                saveas(gcf,'Open_Ellipse.png')
            else
                disp('save images not selected');
            end
            
        else
            disp('Ellipse plot not selected')
        end
        
        %%
        %%%%  Stream Slice %%%%
        if PlotStream == 1
            figure;
            A5=myvec2mat(A4,row,col);
            B5=myvec2mat(B4,row,col);
            u5=myvec2mat(u1,row,col);
            v5=myvec2mat(v1,row,col);
            streamslice(A5,B5,u5,-v5,'noarrows'); %negative v as y axis is flipped to match SALS scan path
            
            set(gca,'Ydir','reverse')
            axis equal
            
            if SavePlots == 1
                
                saveas(gcf,'StreamSlice.fig')
                saveas(gcf,'StreamSlice.png')
            else
                disp('save images not selected');
            end
            
            
        else
            disp('StreamSlice plot not selected')
        end
        
        AvgEccentricity1=mean(ratio);
        %%
        if PlotHistogram == 1
            %% Histograms %%
            %         figure;
            %         edges = [-90:5:90];
            %         h = histogram(Fibreangle,edges);
            %         xlim([-90 90])
            %         xlabel({'Azimuth Angle (?)'});
            %         ylabel({'Frequency'});
            %         saveas(gcf,'hist.png')
            %         saveas(gcf,'hist.fig')
            
            %normalised histogram
            figure;
            edges = (-90:5:90);
            h = histogram(Fibreangle,edges,'Normalization','probability');
            xlim([-90 90])
            ylim([0 1])
            xlabel({'Azimuth Angle (?)'});
            ylabel({'Normalised frequency'});
            
            if SavePlots == 1
                
                saveas(gcf,'hist_norm.png')
                saveas(gcf,'hist_norm.fig')
            else
                disp('save images not selected');
            end
            
        else
            disp('Histogram plot not selected')
        end
        
        %%clear variables to run through for subsequent loops
        clear A2 B2 C2 u2 v2 w2
        
    else
        disp No_Relevant_Files
    end
    
    %     %% Finding the Orientation Index (OI)
    %     %OI = the minimum x distance(angle) containing half the area under the curve
    %     clear x v xq vq1 Area max min
    %     % close all;
    %     clc;
    %
    %     x = 1:1:180; % x values from 1 to 180
    %
    %     % 180 y values from histogram
    %     for j=1:180
    %         v(1,j) = ints(1,j)
    %     end
    %     v = v-min(v);
    %     v = v/max(v);
    %     % v = Hist';
    %     xq=x;
    %
    %
    %     figure;
    %     vq1 = interp1(x,v,xq); % integrate (v,x) at points xq (limits), initially set = x
    %     Area = sum(vq1); % sum to get area under total curve
    %     plot(x,v,'o',xq,vq1,':.');
    %     xlim([0 180]);
    %     title('(Default) Linear Interpolation');
    %
    %     [max, pos] = max(v); %find max peak and its corresponding angle
    %     low = pos-1; %set upper/lower limits based on peak as first guess
    %     up = pos+1;
    %     Area2 = 0; % first guess for OI area = 0
    %
    %     %%
    %     while Area2 < Area/2
    %
    %         xq2 = low:1:up; %integration points to be incremented
    %
    %         vq2 = interp1(x,v,xq2); %integrate smaller region and sum to get area
    %         Area2 = sum(vq2);
    %
    %         low = low-1; %increment limits
    %         up = up+1;
    %
    %     end
    %
    %     figure;
    %     plot(x,v,'o',xq2,vq2,':.');
    %     xlim([0 180]);
    %     title('(Default) Linear Interpolation');
    %     saveas(gcf,'OI.png')
    %     saveas(gcf,'OI.fig')
    %     OI = up-low
    %
    %     Data9 = OI;
    %     xlswrite('Summary.xlsx', Data9, ['I2:I2']);
    %     pause(0.2)
    
    % clear Fibreangle u1 u2 v1 v2 A B A2 B2;
    save SALS_variables.mat
    
end
fprintf(['\nThere are ' num2str(NumberPos) ' positively orientated fibres']);
fprintf(['\nThere are ' num2str(NumberNeg) ' negatively orientated fibres']);
fprintf(['\nThere are ' num2str(NumberZero) ' fibres with an angle of 0 degrees']);