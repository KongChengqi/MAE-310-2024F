clear all; clc; clf; % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

%存储结果
resultL2 = zeros(1,8);
resultH1 = zeros(1,8);
resulth = zeros(1,8);
L2_error=zeros(1,8);
H1_error=zeros(1,8);

%exact solution
exact = @(x)x.^5;%exact solution
exact_square = @(x)x.^10;%exact u square
exact_du = @(x)5*x.^4;
exact_square_du = @(x)25*x.^8;


% Setup the mesh
pp   = 2;              % polynomial degree1
n_en = pp + 1;         % number of element or local nodes
for n_el = 2:2:16     % mesh with element number from 2 to 16
    n_np = n_el * pp + 1;  % number of nodal points
    n_eq = n_np - 1;       % number of equations

    hh = 1.0 / (n_np - 1); % space between two adjacent nodes
    x_coor = 0 : hh : 1;    % number of equations
    
    %set IEN
    IEN = zeros(n_el, n_en);
    for ee = 1 : n_el
        for aa = 1 : n_en
            IEN(ee, aa) = (ee - 1) * pp + aa;
        end
    end
    % Setup the ID array for the problem
    ID = 1 : n_np;
    ID(end) = 0;

    % Setup the quadrature rule
    n_int = 10;
    [xi, weight] = Gauss(n_int, -1, 1);

    % allocate the stiffness matrix
    K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
    F = zeros(n_eq, 1);

    % Assembly of the stiffness matrix and load vector
    for ee = 1 : n_el
        k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
        f_ele = zeros(n_en, 1);    % allocate a zero element load vector
        x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(aa, ee)
        % quadrature loop
        for qua = 1 : n_int    
            dx_dxi = 0.0;
            x_l = 0.0;
            for aa = 1 : n_en
                x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
                dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
            end
            dxi_dx = 1.0 / dx_dxi;
            for aa = 1 : n_en
                f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
                for bb = 1 : n_en
                    k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
                end
            end
        end
        % Assembly of the matrix and vector based on the ID or LM data
        for aa = 1 : n_en
            P = ID(IEN(ee,aa));
            if(P > 0)
                F(P) = F(P) + f_ele(aa);
                for bb = 1 : n_en
                    Q = ID(IEN(ee,bb));
                    if(Q > 0)
                        K(P, Q) = K(P, Q) + k_ele(aa, bb);
                    else
                        F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
                    end
                end
            end
        end
    end
    % ee = 1 F = NA(0)xh
    F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;

    % Solve Kd = F equation
    d_temp = K \ F;
    disp = [d_temp; g];
    % Postprocessing: visualization
    %plot(x_coor, disp, '--r','LineWidth',3);
    %x_sam = 0 : 0.01 : 1;
    %y_sam = x_sam.^5;
    %hold on;
    %plot(x_sam, y_sam, '-k', 'LineWidth', 3);
    n_sam = 20;
    xi_sam = -1 : (2/n_sam) : 1;
    x_sam = zeros(n_el * n_sam + 1, 1);
    y_sam = x_sam; % store the exact solution value at sampling points
    u_sam = x_sam; % store the numerical solution value at sampling pts
    
    for ee = 1 : n_el
        x_ele = x_coor( IEN(ee, :) );
        u_ele = disp( IEN(ee, :) );
        if ee == n_el
            n_sam_end = n_sam+1;
        else
            n_sam_end = n_sam;
        end
        for ll = 1 : n_sam_end
            x_l = 0.0;
            u_l = 0.0;
            for aa = 1 : n_en
                x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
                u_l = u_l + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
            end
            x_sam( (ee-1)*n_sam + ll ) = x_l;
            u_sam( (ee-1)*n_sam + ll ) = u_l;
            y_sam( (ee-1)*n_sam + ll ) = x_l^5;     
        end
    end

    %计算eL2和eH1
    nqp=10;
    [xi,weight]=Gauss(nqp,-1,1);
    nL2=0;
    nH1=0;
    L2_down=0;
    H1_down=0;
    for ee = 1: n_el
        x_ele=x_coor(IEN(ee,:));
        u_ele=disp(IEN(ee,:));

        for ll = 1 : nqp
    x_l = 0.0; uh = 0.0; dx_dxi = 0.0; uh_xi = 0.0;
    for aa = 1 : n_en
      x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(ll), 0);
      uh     = uh     + u_ele(aa) * PolyShape(pp, aa, xi(ll), 0);
      dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(ll), 1);
      uh_xi  = uh_xi  + u_ele(aa) * PolyShape(pp, aa, xi(ll), 1);
    end
    dxi_dx = 1.0 / dx_dxi;
    nL2 = nL2 + weight(ll) * (uh - exact(x_l))^2 * dx_dxi;
    L2_down = L2_down + weight(ll) * exact(x_l)^2 * dx_dxi;
    nH1 = nH1 + weight(ll) * ( uh_xi * dxi_dx - exact_du(x_l) )^2 * dx_dxi;
    H1_down = H1_down + weight(ll) * exact_du(x_l)^2 * dx_dxi;

        end
    end

    %保存
    nH1=nH1^0.5;
    nL2=nL2^0.5;
    H1_down=H1_down^0.5;
    L2_down=L2_down^0.5;
    H1_error(:,n_el/2) = nH1/H1_down;
    L2_error(:,n_el/2) = nL2/L2_down;
    resultH1(:,n_el/2)=log(H1_error(:,n_el/2));
    resultL2(:,n_el/2)=log(L2_error(:,n_el/2));
    resulth(:,n_el/2)=log(hh);
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

















% EOF