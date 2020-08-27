function [lambda,sigma] = weightcal(coord,model,cc,nugget,ktype,index_v1,index_v2)

nst = size(model,1);
model_rotationmatrix = zeros(3,3,nst);
x=coord;  % coordinates of data + locations to estimate
C = zeros(size(x,1),size(x,1));

for i = 1:nst
  A = setrot(model,i);
  model_rotationmatrix(:,:,i) = A;
  R = model_rotationmatrix(:,:,i);
  h = x*R;
  h = h*h';
  h = sqrt(max(0,-2*h+diag(h)*ones(1,size(x,1))+ones(size(x,1),1)*diag(h)'));
  C (:,:,i)= cova(model(i,1),h);
 % C = C + cc(i)*cova(model(i,1),h,b(i));
  %C(1:1+size(C,1):end,i) = 1;
end


ndata=size(coord,1);


if ktype<1 %Simple cokriging 
covar_11 =zeros(index_v1(1)+index_v1(2),index_v1(1)+index_v1(2),nst); 
covar_12 = zeros(index_v1(1)+index_v1(2),index_v2(1)+index_v2(2),nst);
covar_21 = zeros(index_v2(1)+index_v2(2),index_v1(1)+index_v1(2),nst);
covar_22 = zeros(index_v2(1)+index_v2(2),index_v2(1)+index_v2(2),nst);
covar_101 = zeros(index_v1(1)+index_v1(2),1,nst);
covar_201 = zeros(index_v2(1)+index_v2(2),1,nst);
covar_102 = zeros(index_v1(1)+index_v1(2),1,nst);
covar_202 = zeros(index_v2(1)+index_v2(2),1,nst);



%left hand
for i = 1:nst
covar_11(:,:,i)= cc(i,1)*C(1:index_v1(1)+index_v1(2),1:index_v1(1)+index_v1(2),i);
covar_12 (:,:,i)= cc(i,2)*C(1:index_v1(1)+index_v1(2),index_v1(1)+index_v1(2)+1:ndata-1,i);
covar_21(:,:,i) = cc(i,3)*C(index_v1(1)+index_v1(2)+1:ndata-1,1:index_v1(1)+index_v1(2),i);
covar_22(:,:,i)= cc(i,4)*C(index_v1(1)+index_v1(2)+1:ndata-1,index_v1(1)+index_v1(2)+1:ndata-1,i);
% Right hand
covar_101(:,:,i)= cc(i,1)*C(1:index_v1(1)+index_v1(2),ndata,i);
covar_201(:,:,i)= cc(i,2)*C(index_v1(1)+index_v1(2)+1:ndata-1,ndata,i);
covar_102(:,:,i)= cc(i,3)*C(1:index_v1(1)+index_v1(2),ndata,i);
covar_202(:,:,i)= cc(i,4)*C(index_v1(1)+index_v1(2)+1:ndata-1,ndata,i);
end



covar_11 =sum(covar_11,3)+nugget(1,1);
covar_12 = sum(covar_12,3)+nugget(1,2);
covar_21 = sum(covar_21,3)+nugget(1,3);
covar_22 = sum(covar_22,3)+nugget(1,4);
covar_101 = sum(covar_101,3)+nugget(1,1);
covar_201 = sum(covar_201,3)+nugget(1,3);
covar_102 = sum(covar_102,3)+nugget(1,2);
covar_202 = sum(covar_202,3)+nugget(1,4);



% left and Right Covariance Matrix
left = [covar_11 covar_12; covar_21 covar_22];
right = [covar_101 covar_102;covar_201 covar_202];

% Calculation of weigths and variance 
lambda = left\right;
sigma_1 = covar_11(1,1)-[covar_101;covar_201]'*lambda(:,1);
sigma_2 = covar_22(1,1)-[covar_102;covar_202]'*lambda(:,2);
%sigma_3 = covar_21(1,1);
%sigma_4 = covar_12(1,1);
sigma_1 = covar_11(1,1) - lambda'*right - right'*lambda + lambda'*left*lambda;
sigma_3 = covar_21(1,1) - lambda'*right - right'*lambda + lambda'*left*lambda;
sigma_2 = covar_22(1,1) - lambda'*right - right'*lambda + lambda'*left*lambda;



sigma=[sigma_1(1,1) sigma_3(1,2); sigma_3(2,1) sigma_2(2,2)];


elseif ktype<2 %Multicollocated cokriging 
covar_22 =zeros(size(1:index_v2,2),size(1:index_v2,2),nst); 
covar_21 = zeros(size(1:index_v2,2),size(1:index_v1+1,2),nst);
covar_12 = zeros(size(1:index_v1+1,2),size(1:index_v2,2),nst);
covar_11 = zeros(size(1:index_v1+1,2),size(1:index_v1+1,2),nst);
covar_201 = zeros(size(1:index_v2,2),1,nst);
covar_101 = zeros(size(1:index_v1+1,2),1,nst);
covar_202 = zeros(size(1:index_v2,2),1,nst);
covar_102 = zeros(size(1:index_v1+1,2),1,nst);
%left hand
for i = 1:nst
covar_22(:,:,i)= cc(i,4)*C(1:index_v2,1:index_v2,i);
covar_21 (:,:,i)= cc(i,3)*C(1:index_v2,1:index_v1+1,i);
covar_12(:,:,i) = cc(i,2)*C(1:index_v1+1,1:index_v2,i);
covar_11(:,:,i)= cc(i,1)*C(1:index_v1+1,1:index_v1+1,i);
% Right hand
covar_202(:,:,i)= cc(i,4)*C(1:index_v2,ndata,i);
covar_102(:,:,i)= cc(i,3)*C(1:index_v1+1,ndata,i);
covar_201(:,:,i)= cc(i,2)*C(1:index_v2,ndata,i);
covar_101(:,:,i)= cc(i,1)*C(1:index_v1+1,ndata,i);
end


covar_11 =sum(covar_11,3)+nugget(1,1);
covar_12 = sum(covar_12,3)+nugget(1,2);
covar_21 = sum(covar_21,3)+nugget(1,3);
covar_22 = sum(covar_22,3)+nugget(1,4);
covar_101 = sum(covar_101,3)+nugget(1,1);
covar_201 = sum(covar_201,3)+nugget(1,3);
covar_102 = sum(covar_102,3)+nugget(1,2);
covar_202 = sum(covar_202,3)+nugget(1,4);

% left and Right Covariance Matrix
left = [covar_22 covar_21; covar_12 covar_11];
right = [covar_202 covar_201;covar_102 covar_101];
%right = [covar_201;covar_101];



% Calculation of weigths and variance 
lambda = left\right;
%lambda=pinv(left)*right;
sigma_1 = covar_22(1,1)-[covar_202;covar_102]'*lambda(:,1);
sigma_2 = covar_11(1,1)-[covar_201;covar_101]'*lambda(:,2);
%sigma_1 = covar_11(1,1) - lambda'*right - right'*lambda + lambda'*left*lambda;
%sigma_2 = covar_22(1,1) - lambda'*right - right'*lambda + lambda'*left*lambda;
sigma=[sigma_1 sigma_2];   
    
    
elseif ktype<3 %Collocated cokriging 
covar_22 =zeros(size(1:ndata-1,2),size(1:ndata-1,2),nst); 
covar_21 = zeros(size(1:ndata-1,2),1,nst);
covar_12 = zeros(1,size(1:ndata-1,2),nst);
covar_11 = zeros(1,1,nst);
covar_202 = zeros(size(1:ndata-1,2),1,nst);
covar_102 = zeros(1,1,nst);
covar_201 = zeros(size(1:ndata-1,2),1,nst);
covar_101 = zeros(1,1,nst);
%left hand
for i = 1:nst
covar_22(:,:,i)= cc(i,4)*C(1:ndata-1,1:ndata-1,i);
covar_21 (:,:,i)= cc(i,3)*C(1:ndata-1,ndata,i);
covar_12(:,:,i) = cc(i,2)*C(ndata,1:ndata-1,i);
covar_11(:,:,i)= cc(i,1)*C(ndata,ndata,i);
% Right hand
covar_202(:,:,i)= cc(i,4)*C(1:ndata-1,ndata,i);
covar_102(:,:,i)= cc(i,2)*C(ndata,ndata,i);
covar_201(:,:,i)= cc(i,3)*C(1:ndata-1,ndata,i);
covar_101(:,:,i)= cc(i,1)*C(ndata,ndata,i);
end
covar_11 =sum(covar_11,3)+nugget(1,1);
covar_12 = sum(covar_12,3)+nugget(1,2);
covar_21 = sum(covar_21,3)+nugget(1,3);
covar_22 = sum(covar_22,3)+nugget(1,4);
covar_101 = sum(covar_101,3)+nugget(1,1);
covar_201 = sum(covar_201,3)+nugget(1,3);
covar_102 = sum(covar_102,3)+nugget(1,2);
covar_202 = sum(covar_202,3)+nugget(1,4);

% left and Right Covariance Matrix
left = [covar_22 covar_21; covar_12 covar_11];
right = [covar_202 covar_201;covar_102 covar_101];

% Calculation of weigths and variance 
lambda = left\right;
sigma_1 = covar_22(1,1)-[covar_202;covar_102]'*lambda(:,1);
sigma_2 = covar_11(1,1)-[covar_201;covar_101]'*lambda(:,2);


sigma=[sigma_1 sigma_2];



elseif ktype<4 %Simple cokriging 
covar_11 =zeros(size(1:ndata-1,2),size(1:ndata-1,2),nst); 
covar_12 = zeros(size(1:ndata-1,2),size(1:ndata-1,2),nst);
covar_21 = zeros(size(1:ndata-1,2),size(1:ndata-1,2),nst);
covar_22 = zeros(size(1:ndata-1,2),size(1:ndata-1,2),nst);
covar_101 = zeros(size(1:ndata-1,2),size(ndata,2),nst);
covar_201 = zeros(size(1:ndata-1,2),size(ndata,2),nst);
covar_102 = zeros(size(1:ndata-1,2),size(ndata,2),nst);
covar_202 = zeros(size(1:ndata-1,2),size(ndata,2),nst);



%left hand
for i = 1:nst
covar_11(:,:,i)= cc(i,1)*C(1:ndata-1,1:ndata-1,i);
covar_12 (:,:,i)= cc(i,2)*C(1:ndata-1,1:ndata-1,i);
covar_21(:,:,i) = cc(i,3)*C(1:ndata-1,1:ndata-1,i);
covar_22(:,:,i)= cc(i,4)*C(1:ndata-1,1:ndata-1,i);
% Right hand
covar_101(:,:,i)= cc(i,1)*C(1:ndata-1,ndata,i);
covar_201(:,:,i)= cc(i,2)*C(1:ndata-1,ndata,i);
covar_102(:,:,i)= cc(i,3)*C(1:ndata-1,ndata,i);
covar_202(:,:,i)= cc(i,4)*C(1:ndata-1,ndata,i);
end



covar_11 =sum(covar_11,3)+nugget(1,1);
covar_12 = sum(covar_12,3)+nugget(1,2);
covar_21 = sum(covar_21,3)+nugget(1,3);
covar_22 = sum(covar_22,3)+nugget(1,4);
covar_101 = sum(covar_101,3)+nugget(1,1);
covar_201 = sum(covar_201,3)+nugget(1,3);
covar_102 = sum(covar_102,3)+nugget(1,2);
covar_202 = sum(covar_202,3)+nugget(1,4);



% left and Right Covariance Matrix
left = [covar_11 covar_12; covar_21 covar_22];
right = [covar_101 covar_102;covar_201 covar_202];

% Calculation of weigths and variance 
lambda_T = left\right;
sigma_1 = covar_11(1,1)-[covar_101;covar_201]'*lambda_T(:,1);
sigma_2 = covar_22(1,1)-[covar_102;covar_202]'*lambda_T(:,2);

sigma=[sigma_1 sigma_2];



elseif ktype<5 %Simple cokriging 
covar_11 =zeros(size(1:index_v1,2),size(1:index_v1,2),nst); 
covar_12 = zeros(size(1:index_v1,2),size(1:index_v2,2),nst);
covar_21 = zeros(size(1:index_v2,2),size(1:index_v1,2),nst);
covar_22 = zeros(size(1:index_v2,2),size(1:index_v2,2),nst);
covar_101 = zeros(size(1:ndata-1,2),size(ndata,2),nst);
covar_201 = zeros(size(1:ndata-1,2),size(ndata,2),nst);
covar_102 = zeros(size(1:ndata-1,2),size(ndata,2),nst);
covar_202 = zeros(size(1:ndata-1,2),size(ndata,2),nst);



%left hand
for i = 1:nst
covar_11(:,:,i)= cc(i,1)*C(1:ndata-1,1:ndata-1,i);
covar_12 (:,:,i)= cc(i,2)*C(1:ndata-1,1:ndata-1,i);
covar_21(:,:,i) = cc(i,3)*C(1:ndata-1,1:ndata-1,i);
covar_22(:,:,i)= cc(i,4)*C(1:ndata-1,1:ndata-1,i);
% Right hand
covar_101(:,:,i)= cc(i,1)*C(1:ndata-1,ndata,i);
covar_201(:,:,i)= cc(i,2)*C(1:ndata-1,ndata,i);
covar_102(:,:,i)= cc(i,3)*C(1:ndata-1,ndata,i);
covar_202(:,:,i)= cc(i,4)*C(1:ndata-1,ndata,i);
end



covar_11 =sum(covar_11,3)+nugget(1,1);
covar_12 = sum(covar_12,3)+nugget(1,2);
covar_21 = sum(covar_21,3)+nugget(1,3);
covar_22 = sum(covar_22,3)+nugget(1,4);
covar_101 = sum(covar_101,3)+nugget(1,1);
covar_201 = sum(covar_201,3)+nugget(1,3);
covar_102 = sum(covar_102,3)+nugget(1,2);
covar_202 = sum(covar_202,3)+nugget(1,4);



% left and Right Covariance Matrix
left = [covar_11 covar_12; covar_21 covar_22];
right = [covar_101 covar_102;covar_201 covar_202];

% Calculation of weigths and variance 
lambda = left\right;
sigma_1 = covar_11(1,1)-[covar_101;covar_201]'*lambda(:,1);
sigma_2 = covar_22(1,1)-[covar_102;covar_202]'*lambda(:,2);

sigma=[sigma_1 sigma_2];





elseif ktype<6 %Multicollocated cokriging 
covar_22 =zeros(size(1:index_v2(1)+index_v2(2),2),size(1:index_v2(1)+index_v2(2),2),nst); 
covar_21 = zeros(size(1:index_v2(1)+index_v2(2),2),size(1:index_v1(1)+index_v1(2)+1,2),nst);
covar_12 = zeros(size(1:index_v1(1)+index_v1(2)+1,2),size(1:index_v2(1)+index_v2(2),2),nst);
covar_11 = zeros(size(1:index_v1(1)+index_v1(2)+1,2),size(1:index_v1(1)+index_v1(2)+1,2),nst);
covar_201 = zeros(size(1:index_v2(1)+index_v2(2),2),1,nst);
covar_101 = zeros(size(1:index_v1(1)+index_v1(2)+1,2),1,nst);
covar_202 = zeros(size(1:index_v2(1)+index_v2(2),2),1,nst);
covar_102 = zeros(size(1:index_v1(1)+index_v1(2)+1,2),1,nst);
%left hand
for i = 1:nst
covar_22(:,:,i)= cc(i,4)*C(1:index_v2(1)+index_v2(2),1:index_v2(1)+index_v2(2),i);
covar_21 (:,:,i)= cc(i,3)*C(1:index_v2(1)+index_v2(2),index_v2(1)+index_v2(2)+1:ndata,i);
covar_12(:,:,i) = cc(i,2)*C(index_v2(1)+index_v2(2)+1:ndata,1:index_v2(1)+index_v2(2),i);
covar_11(:,:,i)= cc(i,1)*C(index_v2(1)+index_v2(2)+1:ndata,index_v2(1)+index_v2(2)+1:ndata,i);
% Right hand
covar_202(:,:,i)= cc(i,4)*C(1:index_v2(1)+index_v2(2),ndata,i);
covar_102(:,:,i)= cc(i,3)*C(index_v2(1)+index_v2(2)+1:ndata,ndata,i);
covar_201(:,:,i)= cc(i,2)*C(1:index_v2(1)+index_v2(2),ndata,i);
covar_101(:,:,i)= cc(i,1)*C(index_v2(1)+index_v2(2)+1:ndata,ndata,i);
end


covar_11 =sum(covar_11,3)+nugget(1,1);
covar_12 = sum(covar_12,3)+nugget(1,2);
covar_21 = sum(covar_21,3)+nugget(1,3);
covar_22 = sum(covar_22,3)+nugget(1,4);
covar_101 = sum(covar_101,3)+nugget(1,1);
covar_201 = sum(covar_201,3)+nugget(1,3);
covar_102 = sum(covar_102,3)+nugget(1,2);
covar_202 = sum(covar_202,3)+nugget(1,4);

% left and Right Covariance Matrix
left = [covar_22 covar_21; covar_12 covar_11];
right = [covar_202 covar_201;covar_102 covar_101];
%right = [covar_201;covar_101];



% Calculation of weigths and variance 
lambda = left\right;
%lambda=pinv(left)*right;
sigma_1 = covar_22(1,1)-[covar_202;covar_102]'*lambda(:,1);
sigma_2 = covar_11(1,1)-[covar_201;covar_101]'*lambda(:,2);
%sigma_1 = covar_11(1,1) - lambda'*right - right'*lambda + lambda'*left*lambda;
%sigma_2 = covar_22(1,1) - lambda'*right - right'*lambda + lambda'*left*lambda;
sigma=[sigma_1 sigma_2];   






end

end








