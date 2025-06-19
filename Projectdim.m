%define all known values
r_wheel=70;
r_wheel_gear=(41)/(2*pi);
r_crank_gear=(17)/2;
l_3=45;l_4=40;l_2=17.5;
L_1=[62.4 64 64.5 65 69.8 70 70.5 71.5 75.7 76 76.3 76.5];
T_1=[];
T_4=[];
ALPHA4=[];
i=1;
t2=0;t3=pi/2;
%finding values for each ground link val
% ue
while i<=length(L_1)
    % Define tolerance for convergence 
    tol = 0.0001;
    % Assume values for unknowns i.e. t1 and t4
    t1=pi/3;t4=pi/6;
    % Initialize variables for iteration
    dt1=1;dt3=1;   
    count=0;
    % Iteratively solve for t3 and t4 until convergence using Newton rapson method  
    while abs(dt1)>tol || abs(dt4)>tol
        % Define the system of equations
        e1=l_3*cos(t3)+l_4*cos(t4)-l_2*cos(t2)-L_1(i)*cos(t1); % Equation 1
        e2=l_3*sin(t3)+l_4*sin(t4)-l_2*sin(t2)-L_1(i)*sin(t1); % Equation 2
        % Calculate the partial derivatives of e1 and e2 with respect to t1 and t4
        de1_dt1=L_1(i)*sin(t1); % Partial derivative of e1 with respect to t1
        de1_dt4=-l_4*sin(t4); % Partial derivative of e1 with respect to t4
        de2_dt1=-L_1(i)*cos(t1); % Partial derivative of e2 with respect to t1
        de2_dt4=l_4*cos(t4); % Partial derivative of e2 with respect to t4
        % Create the Jacobian matrix
        J=[de1_dt1, de1_dt4; de2_dt1, de2_dt4];
        B=[-e1;-e2];
        X=J\B;
        dt1=X(1);
        dt4=X(2);
        % Update angles
        t1=t1+dt1;
        t4=t4+dt4;
        count=count+1;
        % Break the loop if it doesn't converge
        if count>50
            continue;
        end
    end 
    T_1(end+1)=t1;
    T_4(end+1)=t4;
    i=i+1;
end
%finding velocity values
i=1;
w2=((15*(500/18))*r_wheel_gear)/(r_wheel*r_crank_gear);
%The three arrays which will contain the values of length, angular velocity and angle 
%are defined
L_1g=[];
ALPHA4g=[];
T_1g=[];
%kinematic analysis
while i<=length(L_1)
    Q=[l_3*sin(t3) l_4*sin(T_4(i));l_3*cos(t3) l_4*cos(T_4(i))]\[l_2*sin(t2);l_2*cos(t2)];
    R=[l_3*sin(t3) l_4*sin(T_4(i));l_3*cos(t3) l_4*cos(T_4(i))]\[l_2*cos(t2)-l_3*cos(t3)*Q(1)^2-l_4*cos(T_4(i))*Q(2)^2;-l_2*sin(t2)+l_3*sin(t3)*Q(1)^2+l_4*sin(T_4(i))*Q(2)^2];
    ALPHA4(end+1)=R(2)*(w2^2);
    i=i+1;
end
i=1;
%find values for which the links satisfy Grashof's condition and update
%arrays
while i<=length(L_1)
    if L_1(i)+l_2<l_3+l_4
        L_1g(end+1)=L_1(i);ALPHA4g(end+1)=ALPHA4(i);T_1g(end+1)=T_1(i);
    end
    i=i+1;
end
alphadiff=[];
i=1;
while i<=length(L_1g)
    alphadiff(end+1)=abs(ALPHA4(i)-3);
    i=i+1;
end
i=find(alphadiff==min(alphadiff));
%find the final value of ground link length and angle
l_1=L_1g(i);t1=T_1(i);