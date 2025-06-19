% Define tolerance for convergence 
tol = 0.0001;
i=0;
%start the loop (10 cycles)
while i<=3600   
    % Open
    t2=i*(pi/180);t3=pi/1.5;t4=pi/3;
    % Initialize variables for iteration
    dt1=1;dt3=1;   
    count=0;
    % Iteratively solve for t3 and t4 until convergence using Newton rapson method  
    while abs(dt3)>tol || abs(dt4)>tol
        % Define the system of equations
        e1=l_3*cos(t3)+l_4*cos(t4)-l_2*cos(t2)-l_1*cos(t1); % Equation 1
        e2=l_3*sin(t3)+l_4*sin(t4)-l_2*sin(t2)-l_1*sin(t1); % Equation 2
        % Calculate the partial derivatives of e1 and e2 with respect to t3 and t4
        de1_dt3=-l_3*sin(t3); % Partial derivative of e1 with respect to t3
        de1_dt4=-l_4*sin(t4); % Partial derivative of e1 with respect to t4
        de2_dt3=l_3*cos(t3); % Partial derivative of e2 with respect to t3
        de2_dt4=l_4*cos(t4); % Partial derivative of e2 with respect to t4
        % Create the Jacobian matrix
        J=[de1_dt3, de1_dt4; de2_dt3, de2_dt4];
        B=[-e1;-e2];
        X=J\B;
        dt3=X(1);
        dt4=X(2);
        % Update angles
        t3=t3+dt3;
        t4=t4+dt4;
        count=count+1;
        % Break the loop if it doesn't converge
        if count>50
            continue;
        end
    end 
    %aniimation
    plot([0,l_1*cos(t1)],[0,l_1*sin(t1)], ...
        [0,-l_2*cos(t2)],[0,-l_2*sin(t2)], ...
        [-l_2*cos(t2),-l_2*cos(t2)+l_3*cos(t3)],[-l_2*sin(t2),-l_2*sin(t2)+l_3*sin(t3)], ...
        [l_1*cos(t1),-l_2*cos(t2)+l_3*cos(t3)],[l_1*sin(t1),-l_2*sin(t2)+l_3*sin(t3)], ...
        LineWidth=8);
    axis([-75 75 -75 75]);
    axis square;
    %the simulation runs at a speed that is a factor of the actual angular velocity
    pause((2*pi)/(1000000*w2));
    %clear the graph
    clf;
    i=i+1;
end