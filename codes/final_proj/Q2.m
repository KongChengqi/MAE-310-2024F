clear all; clc;

E=1e9;%杨氏模量
nu=0.3;%泊松比
R=0.5;%孔的尺寸
LL=2;%板的尺寸
Tx=1e4;%受到拉力
fd=2;%自由度
n_en=4;%四边形网格
D=zeros(3,3);
D(1,1)=E/(1-nu^2);D(1,2)=E*nu/(1-nu^2);D(1,3)=0;
D(2,1)=E*nu/(1-nu^2);D(2,2)=E/(1-nu^2);D(2,3)=0;
D(3,1)=0;D(3,2)=0;D(3,3)=E/(2*(1+nu));%建立D矩阵


% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);%先采用了老师的高斯积分

%存储结果
resultL2 = zeros(1,8);
resultH1 = zeros(1,8);
resulth = zeros(1,8);
L2_error=zeros(1,8);
H1_error=zeros(1,8);

for nn=1:1:4
    if nn==1
       [ID, IEN, numNodes] = Getmsh('mesh_1.msh');   

n_eq = numNodes*2;

%LM = ID(IEN);

% allocate the stiffness matrix and load vector
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  k_ele = zeros(fd*n_en, fd*n_en); % element stiffness matrix
  f_ele = zeros(fd*n_en, 1);    % element load vector
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    f(1)=f_x(x_l,y_l);
    f(2)=f_y(x_l,y_l);

    for aa = 1 : n_en
      Na = Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
      B1 = [Na_x , 0;
            0    , Na_y;
             Na_y , Na_x];%写出来B1
  
      for i = 1:fd
      
                f_ele(fd* (aa - 1) + i) = f_ele(fd* (aa - 1) + i) + weight(ll) * detJ * f(i)* Na;
      for bb = 1 : n_en
        Nb = Quad(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        B2 = [Nb_x , 0;
             0    , Nb_y;
             Nb_y , Nb_x];%写出来B2
        
        for j=1:fd
            switch i
                    case 1
                    ei = [1, 0];
                    case 2
                     ei = [0, 1];
            end
            switch j
                    case 1
                    ej = [1;
                          0];
                    case 2
                     ej = [0;
                           1];
             end
             k_ele(fd* (aa - 1) + i, fd* (bb - 1) + j) = k_ele(fd* (aa - 1) + i,fd * (bb - 1) + j) + weight(ll) * detJ * ei* B1' * D * B2 *ej;
        end 
      end
       end % end of bb loop
    end% end of aa loop
  end% end of quadrature loop

  for aa = 1 : n_en
      for i = 1:fd
        PP = ID(IEN(ee,aa),i);
            if PP > 0
              F(PP) = F(PP) + f_ele(fd * (aa - 1) + i);
                for bb = 1 : n_en
                    for j=1:fd
                    QQ = ID(IEN(ee,bb),j);
                    if QQ > 0
                     K(PP, QQ) = K(PP, QQ) + k_ele(fd*(aa - 1) + i, fd*(bb - 1) + j);
                 else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
                    end
                    end
                end
            end  
      end
  end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, fd);

for ii = 1 : n_np
    for j=1:fd
        index = ID(ii,j);
            if index > 0
                disp(ii,j) = dn(index);
            else
                % modify disp with the g data. Here it does nothing because g is zero
            end
    end
end

% save the solution vector and number of elements to disp with name
% disp.mat
save("disp", "disp", "n_el_x", "n_el_y");


%计算error
    nL2=zeros(1,fd);
    nH1=zeros(1,fd);

    %写出精确解
    x=ID(2,:);
    y=ID(2,:);
    r=sqrt(x^2+y^2);
    tan(thetalim)=y/x;
    rr=Tx*(1-R^2/r^2)+Tx/2*(1-4*R^2/r^2+3*R^4/r^4)*cos(2*thetalim);
    thth=Tx/2*(1+R^2/r^2)-Tx/2*(1+3*R^4/r^4)*cos(2*thetalim);
    rth=-Tx/2*(1+2*R^2/r^2-3*R^4/r^4);
   
for ee = 1 : n_el
   for i=1:fd
    x_ele = x_coor(IEN(ee, 1:n_en));
    y_ele = y_coor(IEN(ee, 1:n_en));

    for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    uh = 0.0; uh_xl=0.0;uh_yl=0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
      uh = uh + disp(IEN(ee,aa),i)* Quad(aa, xi(ll), eta(ll));
    end
    exact(1)=u_exact(x_l,y_l);
    exact(2)=v_exact(x_l,y_l);
    exact_x(1)=u_x(x_l,y_l);
    exact_x(2)=v_x(x_l,y_l);
    exact_y(1)=u_y(x_l,y_l);
    exact_y(2)=v_y(x_l,y_l);%终于找到了！这里不小心写串了啊啊
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    for aa = 1 : n_en
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            uh_xl = uh_xl + disp(IEN(ee,aa),i) * (Na_xi * dy_deta - Na_eta * dy_dxi) /detJ;
            uh_yl = uh_yl + disp(IEN(ee,aa),i) * (Na_eta * dx_dxi - Na_xi * dx_deta) /detJ;
    end

        nL2(:,i)= nL2(:,i) + weight(ll) * detJ * (uh - exact(:,i))^2;
        
        nH1(:,i) = nH1(:,i) + weight(ll) * detJ * (( uh_xl - exact_x(:,i) )^2 + ( uh_yl - exact_y(:,i) )^2 );
        
    end
   end
        L2=nL2(:,1)+nL2(:,2);
        H1=nH1(:,1)+nH1(:,2);

end
    H1=H1^0.5;
    L2=L2^0.5;
    H1_error(:,nn) = H1;
    L2_error(:,nn) = L2;
    resultH1(:,nn)=log(H1_error(:,nn));
    resultL2(:,nn)=log(L2_error(:,nn));
    resulth(:,nn)=log(hx);


 else nn==2
 [ID, IEN, numNodes] = Getmsh('mesh_2.msh');  
 n_eq = numNodes*2;

%LM = ID(IEN);

% allocate the stiffness matrix and load vector
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  k_ele = zeros(fd*n_en, fd*n_en); % element stiffness matrix
  f_ele = zeros(fd*n_en, 1);    % element load vector
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    f(1)=f_x(x_l,y_l);
    f(2)=f_y(x_l,y_l);

    for aa = 1 : n_en
      Na = Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
      B1 = [Na_x , 0;
            0    , Na_y;
             Na_y , Na_x];%写出来B1
  
      for i = 1:fd
      
                f_ele(fd* (aa - 1) + i) = f_ele(fd* (aa - 1) + i) + weight(ll) * detJ * f(i)* Na;
      for bb = 1 : n_en
        Nb = Quad(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        B2 = [Nb_x , 0;
             0    , Nb_y;
             Nb_y , Nb_x];%写出来B2
        
        for j=1:fd
            switch i
                    case 1
                    ei = [1, 0];
                    case 2
                     ei = [0, 1];
            end
            switch j
                    case 1
                    ej = [1;
                          0];
                    case 2
                     ej = [0;
                           1];
             end
             k_ele(fd* (aa - 1) + i, fd* (bb - 1) + j) = k_ele(fd* (aa - 1) + i,fd * (bb - 1) + j) + weight(ll) * detJ * ei* B1' * D * B2 *ej;
        end 
      end
       end % end of bb loop
    end% end of aa loop
  end% end of quadrature loop

  for aa = 1 : n_en
      for i = 1:fd
        PP = ID(IEN(ee,aa),i);
            if PP > 0
              F(PP) = F(PP) + f_ele(fd * (aa - 1) + i);
                for bb = 1 : n_en
                    for j=1:fd
                    QQ = ID(IEN(ee,bb),j);
                    if QQ > 0
                     K(PP, QQ) = K(PP, QQ) + k_ele(fd*(aa - 1) + i, fd*(bb - 1) + j);
                 else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
                    end
                    end
                end
            end  
      end
  end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, fd);

for ii = 1 : n_np
    for j=1:fd
        index = ID(ii,j);
            if index > 0
                disp(ii,j) = dn(index);
            else
                % modify disp with the g data. Here it does nothing because g is zero
            end
    end
end

% save the solution vector and number of elements to disp with name
% disp.mat
save("disp", "disp", "n_el_x", "n_el_y");


%计算error
    nL2=zeros(1,fd);
    nH1=zeros(1,fd);

    %写出精确解
    x=ID(2,:);
    y=ID(2,:);
    r=sqrt(x^2+y^2);
    tan(thetalim)=y/x;
    rr=Tx*(1-R^2/r^2)+Tx/2*(1-4*R^2/r^2+3*R^4/r^4)*cos(2*thetalim);
    thth=Tx/2*(1+R^2/r^2)-Tx/2*(1+3*R^4/r^4)*cos(2*thetalim);
    rth=-Tx/2*(1+2*R^2/r^2-3*R^4/r^4);
   
for ee = 1 : n_el
   for i=1:fd
    x_ele = x_coor(IEN(ee, 1:n_en));
    y_ele = y_coor(IEN(ee, 1:n_en));

    for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    uh = 0.0; uh_xl=0.0;uh_yl=0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
      uh = uh + disp(IEN(ee,aa),i)* Quad(aa, xi(ll), eta(ll));
    end
    exact(1)=u_exact(x_l,y_l);
    exact(2)=v_exact(x_l,y_l);
    exact_x(1)=u_x(x_l,y_l);
    exact_x(2)=v_x(x_l,y_l);
    exact_y(1)=u_y(x_l,y_l);
    exact_y(2)=v_y(x_l,y_l);%终于找到了！这里不小心写串了啊啊
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    for aa = 1 : n_en
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            uh_xl = uh_xl + disp(IEN(ee,aa),i) * (Na_xi * dy_deta - Na_eta * dy_dxi) /detJ;
            uh_yl = uh_yl + disp(IEN(ee,aa),i) * (Na_eta * dx_dxi - Na_xi * dx_deta) /detJ;
    end

        nL2(:,i)= nL2(:,i) + weight(ll) * detJ * (uh - exact(:,i))^2;
        
        nH1(:,i) = nH1(:,i) + weight(ll) * detJ * (( uh_xl - exact_x(:,i) )^2 + ( uh_yl - exact_y(:,i) )^2 );
        
    end
   end
        L2=nL2(:,1)+nL2(:,2);
        H1=nH1(:,1)+nH1(:,2);

end
    H1=H1^0.5;
    L2=L2^0.5;
    H1_error(:,nn) = H1;
    L2_error(:,nn) = L2;
    resultH1(:,nn)=log(H1_error(:,nn));
    resultL2(:,nn)=log(L2_error(:,nn));
    resulth(:,nn)=log(hx);


    else nn==3
        [ID, IEN, numNodes] = Getmsh('mesh_3.msh'); 
        n_eq = numNodes*2;

%LM = ID(IEN);

% allocate the stiffness matrix and load vector
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  k_ele = zeros(fd*n_en, fd*n_en); % element stiffness matrix
  f_ele = zeros(fd*n_en, 1);    % element load vector
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    f(1)=f_x(x_l,y_l);
    f(2)=f_y(x_l,y_l);

    for aa = 1 : n_en
      Na = Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
      B1 = [Na_x , 0;
            0    , Na_y;
             Na_y , Na_x];%写出来B1
  
      for i = 1:fd
      
                f_ele(fd* (aa - 1) + i) = f_ele(fd* (aa - 1) + i) + weight(ll) * detJ * f(i)* Na;
      for bb = 1 : n_en
        Nb = Quad(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        B2 = [Nb_x , 0;
             0    , Nb_y;
             Nb_y , Nb_x];%写出来B2
        
        for j=1:fd
            switch i
                    case 1
                    ei = [1, 0];
                    case 2
                     ei = [0, 1];
            end
            switch j
                    case 1
                    ej = [1;
                          0];
                    case 2
                     ej = [0;
                           1];
             end
             k_ele(fd* (aa - 1) + i, fd* (bb - 1) + j) = k_ele(fd* (aa - 1) + i,fd * (bb - 1) + j) + weight(ll) * detJ * ei* B1' * D * B2 *ej;
        end 
      end
       end % end of bb loop
    end% end of aa loop
  end% end of quadrature loop

  for aa = 1 : n_en
      for i = 1:fd
        PP = ID(IEN(ee,aa),i);
            if PP > 0
              F(PP) = F(PP) + f_ele(fd * (aa - 1) + i);
                for bb = 1 : n_en
                    for j=1:fd
                    QQ = ID(IEN(ee,bb),j);
                    if QQ > 0
                     K(PP, QQ) = K(PP, QQ) + k_ele(fd*(aa - 1) + i, fd*(bb - 1) + j);
                 else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
                    end
                    end
                end
            end  
      end
  end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, fd);

for ii = 1 : n_np
    for j=1:fd
        index = ID(ii,j);
            if index > 0
                disp(ii,j) = dn(index);
            else
                % modify disp with the g data. Here it does nothing because g is zero
            end
    end
end

% save the solution vector and number of elements to disp with name
% disp.mat
save("disp", "disp", "n_el_x", "n_el_y");


%计算error
    nL2=zeros(1,fd);
    nH1=zeros(1,fd);

    %写出精确解
    x=ID(2,:);
    y=ID(2,:);
    r=sqrt(x^2+y^2);
    tan(thetalim)=y/x;
    rr=Tx*(1-R^2/r^2)+Tx/2*(1-4*R^2/r^2+3*R^4/r^4)*cos(2*thetalim);
    thth=Tx/2*(1+R^2/r^2)-Tx/2*(1+3*R^4/r^4)*cos(2*thetalim);
    rth=-Tx/2*(1+2*R^2/r^2-3*R^4/r^4);
   
for ee = 1 : n_el
   for i=1:fd
    x_ele = x_coor(IEN(ee, 1:n_en));
    y_ele = y_coor(IEN(ee, 1:n_en));

    for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    uh = 0.0; uh_xl=0.0;uh_yl=0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
      uh = uh + disp(IEN(ee,aa),i)* Quad(aa, xi(ll), eta(ll));
    end
    exact(1)=u_exact(x_l,y_l);
    exact(2)=v_exact(x_l,y_l);
    exact_x(1)=u_x(x_l,y_l);
    exact_x(2)=v_x(x_l,y_l);
    exact_y(1)=u_y(x_l,y_l);
    exact_y(2)=v_y(x_l,y_l);%终于找到了！这里不小心写串了啊啊
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    for aa = 1 : n_en
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            uh_xl = uh_xl + disp(IEN(ee,aa),i) * (Na_xi * dy_deta - Na_eta * dy_dxi) /detJ;
            uh_yl = uh_yl + disp(IEN(ee,aa),i) * (Na_eta * dx_dxi - Na_xi * dx_deta) /detJ;
    end

        nL2(:,i)= nL2(:,i) + weight(ll) * detJ * (uh - exact(:,i))^2;
        
        nH1(:,i) = nH1(:,i) + weight(ll) * detJ * (( uh_xl - exact_x(:,i) )^2 + ( uh_yl - exact_y(:,i) )^2 );
        
    end
   end
        L2=nL2(:,1)+nL2(:,2);
        H1=nH1(:,1)+nH1(:,2);

end
    H1=H1^0.5;
    L2=L2^0.5;
    H1_error(:,nn) = H1;
    L2_error(:,nn) = L2;
    resultH1(:,nn)=log(H1_error(:,nn));
    resultL2(:,nn)=log(L2_error(:,nn));
    resulth(:,nn)=log(hx);

    else nn==4
        [ID, IEN, numNodes] = Getmsh('mesh_4.msh'); 
        n_eq = numNodes*2;

%LM = ID(IEN);

% allocate the stiffness matrix and load vector
K = zeros(n_eq, n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : n_el
  x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );
  
  k_ele = zeros(fd*n_en, fd*n_en); % element stiffness matrix
  f_ele = zeros(fd*n_en, 1);    % element load vector
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    f(1)=f_x(x_l,y_l);
    f(2)=f_y(x_l,y_l);

    for aa = 1 : n_en
      Na = Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
      Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
      B1 = [Na_x , 0;
            0    , Na_y;
             Na_y , Na_x];%写出来B1
  
      for i = 1:fd
      
                f_ele(fd* (aa - 1) + i) = f_ele(fd* (aa - 1) + i) + weight(ll) * detJ * f(i)* Na;
      for bb = 1 : n_en
        Nb = Quad(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        B2 = [Nb_x , 0;
             0    , Nb_y;
             Nb_y , Nb_x];%写出来B2
        
        for j=1:fd
            switch i
                    case 1
                    ei = [1, 0];
                    case 2
                     ei = [0, 1];
            end
            switch j
                    case 1
                    ej = [1;
                          0];
                    case 2
                     ej = [0;
                           1];
             end
             k_ele(fd* (aa - 1) + i, fd* (bb - 1) + j) = k_ele(fd* (aa - 1) + i,fd * (bb - 1) + j) + weight(ll) * detJ * ei* B1' * D * B2 *ej;
        end 
      end
       end % end of bb loop
    end% end of aa loop
  end% end of quadrature loop

  for aa = 1 : n_en
      for i = 1:fd
        PP = ID(IEN(ee,aa),i);
            if PP > 0
              F(PP) = F(PP) + f_ele(fd * (aa - 1) + i);
                for bb = 1 : n_en
                    for j=1:fd
                    QQ = ID(IEN(ee,bb),j);
                    if QQ > 0
                     K(PP, QQ) = K(PP, QQ) + k_ele(fd*(aa - 1) + i, fd*(bb - 1) + j);
                 else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
                    end
                    end
                end
            end  
      end
  end
end

% solve the stiffness matrix
dn = K \ F;

% insert dn back into the vector for all nodes
disp = zeros(n_np, fd);

for ii = 1 : n_np
    for j=1:fd
        index = ID(ii,j);
            if index > 0
                disp(ii,j) = dn(index);
            else
                % modify disp with the g data. Here it does nothing because g is zero
            end
    end
end

% save the solution vector and number of elements to disp with name
% disp.mat
save("disp", "disp", "n_el_x", "n_el_y");


%计算error
    nL2=zeros(1,fd);
    nH1=zeros(1,fd);
   %写出精确解
    x=ID(2,:);
    y=ID(2,:);
    r=sqrt(x^2+y^2);
    tan(thetalim)=y/x;
    rr=Tx*(1-R^2/r^2)+Tx/2*(1-4*R^2/r^2+3*R^4/r^4)*cos(2*thetalim);
    thth=Tx/2*(1+R^2/r^2)-Tx/2*(1+3*R^4/r^4)*cos(2*thetalim);
    rth=-Tx/2*(1+2*R^2/r^2-3*R^4/r^4);
for ee = 1 : n_el
   for i=1:fd
    x_ele = x_coor(IEN(ee, 1:n_en));
    y_ele = y_coor(IEN(ee, 1:n_en));

    for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    uh = 0.0; uh_xl=0.0;uh_yl=0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
      uh = uh + disp(IEN(ee,aa),i)* Quad(aa, xi(ll), eta(ll));
    end
    exact(1)=u_exact(x_l,y_l);
    exact(2)=v_exact(x_l,y_l);
    exact_x(1)=u_x(x_l,y_l);
    exact_x(2)=v_x(x_l,y_l);
    exact_y(1)=u_y(x_l,y_l);
    exact_y(2)=v_y(x_l,y_l);%终于找到了！这里不小心写串了啊啊
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
    for aa = 1 : n_en
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            uh_xl = uh_xl + disp(IEN(ee,aa),i) * (Na_xi * dy_deta - Na_eta * dy_dxi) /detJ;
            uh_yl = uh_yl + disp(IEN(ee,aa),i) * (Na_eta * dx_dxi - Na_xi * dx_deta) /detJ;
    end

        nL2(:,i)= nL2(:,i) + weight(ll) * detJ * (uh - exact(:,i))^2;
        
        nH1(:,i) = nH1(:,i) + weight(ll) * detJ * (( uh_xl - exact_x(:,i) )^2 + ( uh_yl - exact_y(:,i) )^2 );
        
    end
   end
        L2=nL2(:,1)+nL2(:,2);
        H1=nH1(:,1)+nH1(:,2);

end
    H1=H1^0.5;
    L2=L2^0.5;
    H1_error(:,nn) = H1;
    L2_error(:,nn) = L2;
    resultH1(:,nn)=log(H1_error(:,nn));
    resultL2(:,nn)=log(L2_error(:,nn));
    resulth(:,nn)=log(hx);
end
end


%plot error L2 and H1
plot(resulth,resultL2,'-r','LineWidth',3);
xlabel('log(h)');
ylabel('log(Error L2)')
title('Error L2 vs mesh size');
hold on;

figure
plot(resulth,resultH1,'-r','LineWidth',3);
xlabel('log(h)');
ylabel('log(Error H1)');
title('Error H1 vs. Mesh Size');

slope_e_L2 = (resultL2(8)-resultL2(1))/(resulth(8)-resulth(1));
slope_e_H1 = (resultH1(8)-resultH1(1))/(resulth(8)-resulth(1));

