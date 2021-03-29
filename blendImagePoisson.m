function imret = blendImagePoisson(im1, im2, roi, targetPosition)

% input: im1 (background), im2 (foreground), roi (in im2), targetPosition (in im1)

%% TODO: compute blended image
%imret = 255-im1;

%%
[h1,w1,dim1]=size(im1);
[h2,w2,dim2]=size(im2);
im10=zeros(h1,w1);
im20=zeros(h2,w2);
%% 计算roi和targetPosition的区域和边界
for i=1:h1
    for j=1:w1
        in=inpolygon(i,j,targetPosition(:,2),targetPosition(:,1));
        im10(i,j)=in;
    end
end
for i=1:h1
    j1=find(im10(i,:),1);
    j2=find(im10(i,:),1,'last');
    im10(i,j1)=2; 
    im10(i,j2)=2;
end
for i=1:w1
    j1=find(im10(:,i),1);
    j2=find(im10(:,i),1,'last');
    im10(j1,i)=2; 
    im10(j2,i)=2;
end
for i=1:h2
    for j=1:w2
        in=inpolygon(i,j,roi(:,2),roi(:,1));
        im20(i,j)=in;
    end
end
for i=1:h2
    j1=find(im20(i,:),1);
    j2=find(im20(i,:),1,'last');
    im20(i,j1)=2; 
    im20(i,j2)=2;
end
for i=1:w2
    j1=find(im20(:,i),1);
    j2=find(im20(:,i),1,'last');
    im20(j1,i)=2; 
    im20(j2,i)=2;
end
%% 
v1=targetPosition(1,1)-roi(1,1);%X指标
v2=targetPosition(1,2)-roi(1,2);%Y指标
im3=zeros(h1,w1,dim1);
im3=im1;
 for i=1:h2
     for j=1:w2
         if(im20(i,j)~=0)
             im3(i+v2,j+v1,:)=im2(i,j,:);
         end
     end
end
%%  构造AX=B X=L^T*L 求解线性方程组
n=sum(sum(im10~=0));
X=zeros(n,3);% RBG三通道,在边界时也包含于X
XX=zeros(n,3);%对应X中元素在原来矩阵的（b,i，j) i行j列 b是边界情况
A=zeros(n);L=zeros(n);
B=zeros(n,3);
m=0;
for i=1:h2
    for j=1:w2
        if(im20(i,j)~=0)
            m=m+1;
            XX(m,2)=i;
            XX(m,3)=j;
            XX(m,1)=im20(XX(m,2),XX(m,3));
        end
    end
end
[FX,FY,FZ]=gradient(double(im3));
[FX_X,FX_Y,FX_Z]=gradient(double(FX));
[FY_X,FY_Y,FY_Z]=gradient(double(FY));
divI=FX_X+FY_Y;
for i=1:n
    if(XX(i,1)==2)
        A(i,i)=1; 
        B(i,:)=im1(XX(i,2)+v2,XX(i,3)+v1,:);
        %continue;
    else
        A(i,i)=-4;
        i0=XX(i,2);j0=XX(i,3);
        j1=find(XX(:,2)==i0-1&XX(:,3)==j0,1);
        A(i,j1)=1;
        j2=find(XX(:,2)==i0+1&XX(:,3)==j0,1);
        A(i,j2)=1;
        j3=find(XX(:,2)==i0&XX(:,3)==j0-1,1);
        A(i,j3)=1;
        j4=find(XX(:,2)==i0&XX(:,3)==j0+1,1);
        A(i,j4)=1;
        B(i,:)=divI(i0+v2,j0+v1,:);
    end
end
X=A\B;
imret=im1;
XX(:,2)=XX(:,2)+v2;
XX(:,3)=XX(:,3)+v1;
for i=1:n
    i0=XX(i,2);
    j0=XX(i,3);
    imret(i0,j0,:)=X(i,:);
end
%imret=im3;








