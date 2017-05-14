function Distance=Distance2(Point1,Point2)
Distance=sqrt((Point1(:,1)-Point2(:,1)).^2+(Point1(:,2)-Point2(:,2)).^2);