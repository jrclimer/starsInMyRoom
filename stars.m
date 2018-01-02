% A = asimuth
% a = altitude
% alpha – right ascension
% delta – declination
% h – hour angle
% lambdao – observer's longitude
% phio – observer's latitude
% thetaL – local sidereal time
% thetaG – Greenwich sidereal time
latitude = 42.3601;% N
longitude = -71.0578;% W
localSiderealTime = 13+40/60+42.1262/60/60;
hourAngle = (localSiderealTime-hygdatav3.ra)/24*360;
A = mod(atan2d(...
    sind(hourAngle)...
    ,cosd(hourAngle)*sind(latitude)-tand(hygdatav3.dec).*cosd(latitude)...
    ),360);
a = asind(sind(latitude)*sind(hygdatav3.dec)+cosd(latitude)*cosd(hygdatav3.dec).*cosd(hourAngle));
goods = hygdatav3.mag<quantile(hygdatav3.mag(a>0),200/sum(a>0))&a>0;
bigDipper = ismember(hygdatav3.proper,{'Alkaid','Mizar','Alioth','Megrez','Phad','Dubhe','Merak'});
polaris = ismember(hygdatav3.proper,'Polaris');
A = -A;

% In A, 0 is south and corresponds with wall III. 90
% A = 0:15:360;
% A = A(:);
% a = ones(size(A));
% A = [ones(15,1)*-135;ones(15,1)*45];
% a = [linspace(0,90,15)';linspace(90,0,15)'];

% A = 180;
% a = 1;
% goods = true(size(A));

width = 9*12+6.5;
height = 24;
length = 12*12+11;
orientation = 90-60.04;% If 0, A = 0 lines up with the center of wall 3. Note that angles
% proceed such that A = -90 = east
orientation = orientation-90;

pos = [...
    -width/2 length/2 0;...
    width/2 length/2 0;...
    width/2 -length/2 0;...
    -width/2 -length/2 0;...
    -width/2 -length/2 height;...
    -width/2 length/2 height;...
    width/2 length/2 height;...
    width/2 -length/2 height]*...
    [cosd(orientation) -sind(orientation) 0;...
    sind(orientation) cosd(orientation) 0;...
    0 0 1];

r = max(sqrt(sum(pos.^2,2)));
starPos = r*[sind(90-a).*cosd(A) sind(90-a).*sind(A) cosd(90-a)];

% B - normal vectors of the planes
n = [cross(diff(pos([1 6],:))',diff(pos([1 2],:))')...I
    cross(diff(pos([2 7],:))',diff(pos([2 3],:))')...II
    cross(diff(pos([3 8],:))',diff(pos([3 4],:))')...III
    cross(diff(pos([4 5],:))',diff(pos([4 1],:))')...IV
    cross(diff(pos([8 7],:))',diff(pos([8 5],:))')...V
    ];
n=n./sqrt(sum(n.^2));
V0 = pos([1 3 4 4 5],:)';% Points on each plane
% u - unit vectors along each star position
u = starPos(goods,:)';
u = u./sqrt(sum(u.^2));

intersectionTest = zeros(5,sum(goods));
for i=1:5
   for j=1:sum(goods)
       intersectionTest(i,j) = dot(n(:,i),u(:,j))>0;
   end
end

% Wall SizeOfStar 
intersectionPoint = NaN(sum(goods),3,5);
for i=1:5
    for j=1:sum(goods)
        if intersectionTest(i,j)
            intersectionPoint(j,:,i) = ...
                dot(n(:,i),V0(:,i))/dot(n(:,i),u(:,j))*u(:,j);
        end
    end
end

intersectionPoint0 = intersectionPoint;

% Now, shift to relative to upper left (SE for roof) corner of plane
for i=1:5
    % First translate, then rotate about corner
    intersectionPoint(:,:,i)=(intersectionPoint(:,:,i)-V0(:,i)')...
        *[cosd(-orientation) -sind(-orientation) 0;...
    sind(-orientation) cosd(-orientation) 0;...
    0 0 1];
end

% Remove stars beyond extent of each wall
intersectionPoint(intersectionPoint(:,1,1)>width|intersectionPoint(:,1,1)<0|intersectionPoint(:,3,1)>height,:,1)=NaN;
intersectionPoint(intersectionPoint(:,2,2)>length|intersectionPoint(:,2,2)<0|intersectionPoint(:,3,2)>height,:,2)=NaN;
intersectionPoint(intersectionPoint(:,1,3)>width|intersectionPoint(:,1,3)<0|intersectionPoint(:,3,3)>height,:,3)=NaN;
intersectionPoint(intersectionPoint(:,2,4)>length|intersectionPoint(:,2,4)<0|intersectionPoint(:,3,4)>height,:,4)=NaN;
intersectionPoint(intersectionPoint(:,1,5)>width|intersectionPoint(:,1,5)<0|intersectionPoint(:,2,5)>length|intersectionPoint(:,2,5)<0,:,5)=NaN;

goods = find(goods);

sizeThresh = quantile(hygdatav3.mag(goods),[25 100]/200);

% Wall StarSize(1,2,3) x-distance y-distance Proper
out = cell(numel(goods),5);
for i=1:numel(goods)
    out{i,1} = find(~isnan(intersectionPoint(i,1,:)),1);
    switch out{i,1}    
        case 3
            out{i,2} = width - intersectionPoint(i,1,out{i,1});
        case {1,5}
                out{i,2} = intersectionPoint(i,1,out{i,1});
        case 2
            out{i,2} = length - intersectionPoint(i,2,out{i,1});
        otherwise
            out{i,2} = intersectionPoint(i,2,out{i,1});
    end
    switch out{i,1}
        case {1,2,3,4}
            out{i,3} = intersectionPoint(i,3,out{i,1});
        case 5
            out{i,3} = intersectionPoint(i,2,out{i,1});
    end
    out{i,4} = sum(hygdatav3.mag(goods(i))>sizeThresh)+1;
    out{i,5} = hygdatav3.proper{goods(i)};
end

out