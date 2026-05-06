\\ Load external function
install("ellheegner_z","GL","ellanal");

\\ Solve linear system of equations for variables XY
solve_xy(S,XY) = {
    my(C,M,n,m,res);
    \\ Extract coefficients of linear terms from S
    M = matrix(2,2,n,m,polcoef(S[n],1,XY[m]));
    \\ Calculate constant term
    C = (M*XY~) - S~;
    \\ Solve the linear system
    res = matsolve(M,C)~;
    return(res);
}

\\ Compute first and second derivatives of S with respect to variables
get_heegner_4d_derivs(S,X,X0) = {
    my(x0,y0,z0,Y1,Y2,Z1,Z2,S1,S2,YZ1,YZ2,S0);
    \\ Define symbolic variables
    [Y1,Y2,Z1,Z2] = ['Y1,'Y2,'Z1,'Z2];
    
    \\ First derivative of S
    S1 = deriv(S,X[2]);
    S1 += deriv(S,X[3])*Y1;
    S1 += deriv(S,X[4])*Z1;
    S0 = substvec(S1,X,X0);
    YZ1 = solve_xy(S0,[Y1,Z1]);
    
    \\ Second derivative of S
    S2 =  deriv(S1,X[2]);
    S2 += deriv(S1,X[3])*Y1;
    S2 += deriv(S1,X[4])*Z1;
    S2 += deriv(S1,Y1)*Y2;
    S2 += deriv(S1,Z1)*Z2;
    S0 = substvec(S2,X,X0);
    S0 = substvec(S0,[Y1,Z1],YZ1);
    YZ2 = solve_xy(S0,[Y2,Z2]);
    
    my(Y,Y0,n,i,j);
    Y = [Y1,Y2,Z1,Z2];
    Y0 = [YZ1,YZ2];

    \\ Extract coefficients from derivatives
    for(n = 1,4, 
        for(i = 1,2, 
            for(j = 1,2, 
                Y0[i][j] = polcoef(Y0[i][j],0,X[n]);
                Y0[i][j] = polcoef(Y0[i][j],0,Y[n]);
            );
        );
    );

    [YZ1,YZ2] = Y0;

    \\ Return derivatives: first order and second order
    return([YZ1[1],YZ1[2],YZ2[1],YZ2[2]]);
}

\\ Get 4x4 transformation matrix for Heegner point computation
get_heegner_M4(S,X,X0,B) = {
    \\ X = [x1,x2,x3,x4]
    \\ X0 = [x0,y0,z0] <- [1,x2,x3,x4]/
    my(x0,y0,z0,S1,y1,y2,z1,z2,e1);
    [x0,y0,z0] = [X0[2],X0[3],X0[4]];
    
    \\ Substitute x1 = 1 into S
    S1 = subst(S,X[1],1);
    [y1,z1,y2,z2] = get_heegner_4d_derivs(S1,X,X0);
    
    my(B2,B3,y10,z10);
    e = z2/y2;
    B2 = B^2; B3 = B^3;
    y10 = (y1*x0-y0);
    z10 = (z1*x0-z0);
    
    \\ Construct the 4x4 transformation matrix
    M4 = [
    1,-x0*B,y10*B2,(-e*y10+z10)*B3;
    0,    B, -y1*B2,   (e*y1-z1)*B3;
    0,    0,    B2,  -e*B3;
    0,    0, 0,     B3
    ];
    return(M4);
}

\\ Get coefficient matrix from a quadric equation
get_mat_4(s1,X) = {
    my(f1,n,m,M1);
    M1 = matrix(4);
    
    \\ Compute second partial derivatives to get the quadratic form matrix
    for(n = 1,4, 
        f1 = deriv(s1,X[n]);
        for(m = 1,4, 
            M1[n,m] = deriv(f1,X[m])/2;
        );
    );
    
    \\ Substitute symbolic variables with 1
    M1 = substvec(M1,X,[1,1,1,1]);
    return(M1);
}

\\ Map a point on a quartic to the one on an elliptic curve
get_tu(arr,P) = {
    my(a,b,c,d,e,x,y,g6,t,u,g4);
    x = P[1]; y = P[2];

    a = arr[1]; b = arr[2]; c = arr[3]; d = arr[4]; e = arr[5];

    \\ Compute degree 4 polynomial g4
    g4 =  (3*b^2-8*a*c)*x^4+4*(b*c-6*a*d)*x^3;
    g4 += 2*(2*c^2-24*a*e-3*b*d)*x^2+4*(c*d-6*b*e)*x+(3*d^2-8*c*e);

    \\ Compute degree 6 polynomial g6
    g6 =  (b^3+8*a^2*d-4*a*b*c)*x^6+2*(16*a^2*e+2*a*b*d-4*a*c^2+b^2*c)*x^5;
    g6 += 5*(8*a*b*e+b^2*d-4*a*c*d)*x^4+20*(b^2*e-a*d^2)*x^3-5*(8*a*d*e+b*d^2-4*b*c*e)*x^2;
    g6 += -2*(16*a*e^2+2*b*d*e-4*c^2*e+c*d^2)*x-(d^3+8*b*e^2-4*c*d*e);
    

    \\ Compute t and u coordinates
    t = 3*g4/4/y^2;
    u = 27*g6/8/y^3;

    return([t,u]);
}

\\ Compute I and J invariants for a quartic
\\ arr: coefficients of a quartic
get_IJ(arr) = {
    my(a,b,c,d,e);
    a = arr[1];
    b = arr[2];
    c = arr[3];
    d = arr[4];
    e = arr[5];
    
    \\ Compute I invariant
    I1 = 12*a*e + (c^2 - 3*d*b);
    \\ Compute J invariant
    J  = (8*a*c - 3*b^2)*9*e + 9*d*(b*c-3*a*d) - 2*c^3;
    
    return([I1,J]);
}

\\ Extract real values from a complex vector, ignoring complex values
get_real_values(Z) = {
    my(Z0,n);
    Z0 = List();
    
    \\ Loop through all elements and keep only real ones
    for(n = 1,#Z, 
        if(abs(imag(Z[n]))<1e-40, 
            listput(Z0,real(Z[n]));
        );
    );
    return(Vec(Z0));
}

\\ Map a point on cubic to a point on quartic
map_cubic_to_quartic(M1,M2,P) = {
    my(S,e1,P2,X);
    \\ Compute determinant of M1*x+M2 to get quartic
    S = matdet(M1*x+M2);
    S = substvec(S,['x1,'x2,'x3,'x4],[1,1,1,1]);
    IJ = get_IJ(Vec(S));
    IJ = substvec(IJ,['x1,'x2,'x3,'x4],[1,1,1,1]);
    
    \\ Initialize elliptic curve from IJ invariants
    e1 = ellinit([-27*IJ[1],-27*IJ[2]]);
    e2 = ellminimalmodel(e1,&tran);
    P2 = ellchangepointinv(P,tran);
    \\print("P2 = ",P2);

    my(t,u,s1,s2,X,x1,y1,n);
    [t,u] = get_tu(Vec(S),[x,y]);
    \\ Solve for x coordinate from the transformed equation
    s1 = t*y^2-S*P2[1];
    s1 = substvec(s1,['x1,'x2,'x3,'x4],[1,1,1,1]);
    s1 = polcoef(s1,0,y);
    X = polroots(s1);
    X = get_real_values(X);
    my(XY);
    XY = List();
    
    \\ For each real root, compute corresponding y coordinate
    for(n = 1,#X, 
        x1 = X[n];
        y1 = sqrt(subst(S,x,x1));
        
        \\ Validate realness of x1 and y1
        if(abs(imag(x1))>1e-40, 
            print("wrong x1!!!!");
            next();
        );
        if(abs(imag(y1))>1e-40, 
            print("wrong y1!!!!");
            next();
        );
        listput(XY,[real(x1),real(y1)]);
    );
    return(XY);
}

\\ Compute F1 and F2 polynomials for mapping to quartic
get_F1F2(M1,M2,P) = {
    my(m1,m2,s,S,X1);
    s = matdet(M1*x+M2);
    \\print("s = ",s);
    S = Vec(s);
    m1 = matadjoint(M1);
    m2 = matadjoint(M2);

    M0 = matadjoint(x*m1+m2);
    
    T = [
    matrix(4),matrix(4),
    matrix(4),matrix(4)
    ];

    \\ Extract coefficients of adjoint matrix
    for(n = 1,4, 
        for(m = 1,4, 
            V = Vec(M0[n,m]);
            for(k = 1,4, 
                T[k][n,m] = polcoef(M0[n,m],4-k,x);
            );
        );
    );
    
    d1 = T[2]/S[1];
    d2 = T[3]/S[5];
    
    my(F1,F2);
    X0 = matrix(4,1,n,m,P[n]);
    X1 = mattranspose(X0);
    F1 = (X1*d1*X0)[1,1];
    F2 = (X1*d2*X0)[1,1];
    return([F1,F2]);
    
}

\\ Map a real-value point (x,y) to (1,x0,y0,z0) in 4-descent
map_quartic_to_4descent(M1,M2,P) = {
    my(X,s1,s2,s3,F1,F2,x0,y0,z0);
    x0 = 'x0; y0 = 'y0; z0 = 'z0;
    x1 = 'x1;

    X = [1,x0,y0,z0];
    \\ Define quadrics s1 and s2
    s1 = X*M1*X~;
    s2 = X*M2*X~;
    X1 = [x1,x0,y0,z0];

    [F1,F2] = get_F1F2(M1,M2,X);
    \\ x1 = -F1/F2
    s3 = F1 + F2*P[1];

    my(s4,s5,s6,A0,A,count,An,an,s5_);
    \\ Resultants to eliminate variables
    s4 = polresultant(s1,s2,x0);
    s5 = polresultant(s1,s3,x0);
    A  = [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15];
    A0 = List();
    count = 0;
    count = 0;
    s5_ = 0;
    
    \\ Build polynomial in z0
    for(i = 0,4, 
        An = polcoef(s5,i,y0);
        an = 0;
        for(j = 0,4-i, 
            count += 1;
            an += A[count]*z0^j;
            listput(A0,polcoef(An,j,z0));
        );
        s5_ += an*y0^i;
    );
    
    \\ Final resultant to get polynomial in z0
    s6 = polresultant(s4,s5_,y0);
    s6 = substvec(s6,A,Vec(A0));
    s6 = substvec(s6,['x1,'x2,'x3,'x4],[1,1,1,1]);
    \\ Computing z0
    
    my(Z,Z0,YZ0,Y,k1,k2,XYZ0);
    Z = polroots(s6);
    Z0 = get_real_values(Z);
    YZ0 = List();
    
    \\ For each real z0 value, find corresponding y0 values
    for(n = 1,#Z0, 
        z0_ = Z0[n];
        S = subst([s4,s5],z0,z0_);
        S = substvec(S,['x1,'x2,'x3,'x4],[1,1,1,1]);
        Y = polroots(S[1]);
        Y = get_real_values(Y);
        for(i = 1,#Y, 
            k1 = abs(subst(S[2],y0,Y[i])); 
            if(k1 < 1e-40, 
                listput(YZ0,[Y[i],z0_]);
            );
        );
    );
    
    XYZ0 = List();
    \\ For each (y0, z0) pair, find corresponding x0 values
    for(n = 1,#YZ0, 
        [y0_,z0_] = YZ0[n];
        S = substvec([s1,s2,s3],[y0,z0],[y0_,z0_]);
        S = substvec(S,['x1,'x2,'x3,'x4],[1,1,1,1]);
        X = polroots(S[1]);
        X = get_real_values(X);
        for(i = 1,#X, 
            k1 = abs(subst(S[2],x0,X[i])); 
            k2 = abs(subst(S[3],x0,X[i])); 
            if(k1 < 1e-40 && k2 < 1e-40, 
                listput(XYZ0,[1,X[i],y0_,z0_]);
            );
        );
    );
    XYZ0 = real(Vec(XYZ0));
    print("there are ",#XYZ0," solutions");
    return(XYZ0);
}

\\ Find Heegner point by z-coordinate using 4-descent
ellheegner_4descent_point_by_z(E,z,M1,M2,prec=400) = {
    \\ E should be a minimalmodel
    my(P,k,XY,y1,S,XYZ0,xn,S1,S2,X,x1_);
    my(M4,Xn);
    P = ellztopoint(E,z);
    \\print("P = ",P);
    
    XY = map_cubic_to_quartic(M1,M2,P);
    my(k,m);
    my(found);
    found = 0;
    
    \\ For each mapped point on quartic, try to find rational point
    for(k = 1,#XY, 
    
        [x1_,y1] = XY[k];
        \\print("quartic: ",[x1_,y1]," there are ",#XY," solutions");
        S = matdet(M1*x+M2);
        
        XYZ0 = map_quartic_to_4descent(M1,M2,[x1_,y1]);
        my(x1,x2,x3,x4);
        [x1,x2,x3,x4] = ['x1,'x2,'x3,'x4];
        X = [x1,x2,x3,x4];
        S1 = X*M1*X~;
        S2 = X*M2*X~;

        my(n1 = floor(prec/40),n);
        
        \\ Try different precision levels
        for(m = 1,#XYZ0, 
            \\print("k, m = ",[k,m]);
            for(n = 1,n1+2, 
                M4 = get_heegner_M4([S1,S2],X,XYZ0[m],10^(40*n));
                Xn = qflll(M4~);
                if(#Xn<4, 
                    next();
                );
                
                \\ Test if we found a rational point
                xn = (Xn*([1,0,0,0]~))~;
                s1 = xn*M1*xn~;
                s2 = xn*M2*xn~;
                if(s1==0&&s2==0, 
                    found = 1;
                    print([k,m,n,xn]);
                    return([found,xn]);
                );
            );
            if(found, 
                break();
            );
        );
    );

    if(!found,
        print("Solution not found!!!");
    );
    
    return([found,-1]);
}

\\ Main function for Heegner point 4-descent algorithm
ellheegner_4descent(E,M1,M2,height,Z=-1) = {
    my(prec1,z1,t1,indx,z2,t2);
    \\ Calculate required precision based on height
    prec1 = floor(height/log(10)/12) + 32;
    print("using precision ",prec1);
    prec2 = floor(prec1*log(10)/log(2));
    print("using bit precision ",prec2);
    
    \\ Compute Heegner point z-coordinate if not provided
    if(Z == -1, 
        localbitprec(prec2);
        Z = ellheegner_z(E,prec2);
    );
    print("Z = ",Z);
    
    localprec(prec1*3);
    z1   = Z[1]; 
    t1   = Z[2][1];
    t2   = Z[2][2];
    indx = round(Z[3]/Z[1]);

    my(m1,found,P,m2);
    found = 0;


    my(disc);
    disc = E.disc;
    print("indx = ",indx);
    
    \\ Try different combinations of m1 and m2
    for(m2 = 0,indx, 
        m1 = indx - m2;
        print("m1 = ",m1);
        z2 = (-z1+m1*t1)/indx;
        res = ellheegner_4descent_point_by_z(E,z2,M1,M2,prec1);
        found = res[1];
        P = res[2];
        if(found, 
            break();
        );

        if(disc>0, 
            print("disc > 0");
            \\ For curves with positive discriminant, try shifted z
            \\ z2 = gadd(z2,gmul2n(Oim,-1))
            z2 += t2/2;
            res = ellheegner_4descent_point_by_z(E,z2,M1,M2,prec1);
            found = res[1];
            P = res[2];
            if(found, 
                break();
            );

        );
    
    );
    /*
    if(!found, 
        for(m1 = 0,indx, 
            z2 = (-z1+m1*t1)/indx;
            if(disc>0, 
                print("disc > 0");
                \\ z2 = gadd(z2,gmul2n(Oim,-1))
                z2 -= t2/2;
                res = ellheegner_4descent_point_by_z(E,z2,M1,M2,prec1);
                found = res[1];
                P = res[2];
                if(found, 
                    break();
                );

            );
        );
    );
    */
    print("P = ",P);
    return(P);
}

\\ Convert a point from original model to minimal model
get_minimalmodel_point(M1,M2,P) = {
    localprec(38);
    my(x1,y1,F,IJ,e1,e2);
    \\ Compute determinant to get quartic
    S  = matdet(x*M1+M2);
    S  = substvec(S,['x1,'x2,'x3,'x4],[1,1,1,1]);
    F  = get_F1F2(M1,M2,P);
    x1 = -F[1]/F[2];
    issquare(subst(S,x,x1),&y1);
    my(P1,P2,tran);
    \\ Compute IJ invariants
    IJ = get_IJ(Vec(S));
    e1 = ellinit([-27*IJ[1],-27*IJ[2]]);
    \\ Transform point
    P1 = get_tu(Vec(S),[x1,y1]);
    e2 = ellminimalmodel(e1,&tran);
    P2 = ellchangepoint(P1,tran);
    P2 = substvec(P2,['x1,'x2,'x3,'x4],[1,1,1,1]);
    print(ellisoncurve(e2,P2));
    print(e2);
    print("z2 = ",ellpointtoz(e2,ellneg(e2,P2)));
    print("height = ",ellheight(e2,P2));
    
    return(P2);
}
