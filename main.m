clc;
clear;
close all;
format compact;

%given variables
a=zeros(377);
const=zeros(377,1);
t_old=zeros(377,1);
t_inf=20;
k_c=50;
k_s=5;
h_inf=500+43/2; %h_inf=500+A/2
qdot=10^7;
dx=0.002;
dy=0.002;
k_average=(k_s+k_c)/2;
t_base=33;
t_limit=85;

for row=1:13
    for col=0:28
        n=row+(col*13);
        if (row==1) %bottom row nodes
           if (col==0) %node 1
               a(n,n+1)=1;
               a(n,n+13)=1;
               a(n,n)=-2;
               const(n,1)=0;
           end
           if (col==28) %node 365
               a(n,n+1)=1;
               a(n,n-13)=1;
               a(n,n)=-2;
               const(n,1)=0;
           end
           if (col>0) && (col<28) %nodes 14:13:352
               if (n==53) || (n==313) %nodes 53&&313 (connected to the base)
                   a(n,n)=1;
                   const(n,1)=t_base;
               else 
                   a(n,n+1)=2;
                   a(n,n+13)=1;
                   a(n,n-13)=1;
                   a(n,n)=-4;
                   const(n,1)=0;
               end
           end
        end
        if ((col==0) || (col==28)) && ((row>1)&&(row<13))
            if (col==0) %node 2:12
                a(n,n+1)=1;
                a(n,n-1)=1;
                a(n,n+13)=2;
                a(n,n)=-4;
                const(n,1)=0;
            else %node 366:376
                a(n,n+1)=1;
                a(n,n-1)=1;
                a(n,n-13)=2;
                a(n,n)=-4;
                const(n,1)=0;
            end
        end
        if (row>1) && (row<13)
           if (row>=9) && ((col>8) && (col<20))
              if col==9 || col==19
                  if row==9
                      if n==126 %node 126
                        a(n,n+1)=k_average;
                        a(n,n-1)=k_s;
                        a(n,n+13)=k_average;
                        a(n,n-13)=k_s;
                        a(n,n)=-2*(k_s+k_average);
                        const(n,1)=-(qdot*dx*dy)/4;
                      else %node 256
                        a(n,n+1)=k_average;
                        a(n,n-1)=k_s;
                        a(n,n-13)=k_average;
                        a(n,n+13)=k_s;
                        a(n,n)=-2*(k_s+k_average);
                        const(n,1)=-(qdot*dx*dy)/4;
                      end
                  else
                      if col==9 %node 127:129
                        a(n,n+1)=k_average;
                        a(n,n-1)=k_average;
                        a(n,n+13)=k_c;
                        a(n,n-13)=k_s;
                        a(n,n)=-(k_s+k_c+2*k_average);
                        const(n,1)=-qdot*(dx*dy)/2;
                      else %257:259
                        a(n,n+1)=k_average;
                        a(n,n-1)=k_average;
                        a(n,n+13)=k_s;
                        a(n,n-13)=k_c;
                        a(n,n)=-(k_s+k_c+2*k_average);
                        const(n,1)=-qdot*(dx*dy)/2;
                      end
                  end
              else
                  if (row == 9) %node 139:13:243
                    a(n,n+1)=k_c;
                    a(n,n-1)=k_s;
                    a(n,n+13)=k_average;
                    a(n,n-13)=k_average;
                    a(n,n)=-(2*k_average+k_c+k_s);
                    const(n,1)=-qdot*dx*dy/2;
                  else %on the chip
                    a(n,n+1)=1;
                    a(n,n-1)=1;
                    a(n,n+13)=1;
                    a(n,n-13)=1;
                    a(n,n)=-4;
                    const(n,1)=-qdot*dx*dy/k_c;
                  end
              end
           else %on the box
               if (col>0) && (col<28)
                  a(n,n+1)=1;
                  a(n,n-1)=1;
                  a(n,n+13)=1;
                  a(n,n-13)=1;
                  a(n,n)=-4;
                  const(n,1)=0;
               end
           end
        end
        if (row==13) %top row nodes
            if (col==0) %node 13
                a(n,n-1)=1;
                a(n,n+13)=1;
                a(n,n)=-(2+(h_inf*dx/k_s));
                const(n,1)=-(h_inf*dx*t_inf/k_s);
            end
            if (col==28) %node 377
                a(n,n-1)=1;
                a(n,n-13)=1;
                a(n,n)=-(2+(h_inf*dx/k_s));
                const(n,1)=-(h_inf*dx*t_inf/k_s);
            end
            if (col>0)&&(col<28)
                if col==9  %node 130
                    a(n,n-1)=2*k_average;
                    a(n,n+13)=k_c;
                    a(n,n-13)=k_s;
                    a(n,n)=-(k_s+k_c+(2*k_average)+(2*h_inf*dx));
                    const(n,1)=-((qdot*dx*dy/2)+(2*h_inf*t_inf*dx));
                end
                if col==19 %node 260
                    a(n,n-1)=2*k_average;
                    a(n,n+13)=k_s;
                    a(n,n-13)=k_c;
                    a(n,n)=-(k_s+k_c+(2*k_average)+(2*h_inf*dx));
                    const(n,1)=-((qdot*dx*dy/2)+(2*h_inf*t_inf*dx));
                end
                if (col>9) && (col<19) %node 143:13:247
                    a(n,n-1)=2;
                    a(n,n+13)=1;
                    a(n,n-13)=1;
                    a(n,n)=-(4+(2*h_inf*dx/k_c));
                    const(n,1)=-((qdot*dx*dy)+(2*h_inf*dx*t_inf))/k_c;
                end
                if (col<9)||(col>19) %node 26:13:117 && node 273:13:364
                    a(n,n-1)=2;
                    a(n,n+13)=1;
                    a(n,n-13)=1;
                    a(n,n)=-(4+(2*h_inf*dx/k_s));
                    const(n,1)=-(2*h_inf*t_inf*dx)/k_s;
                end
            end
        end
    end
end

t_old=a\const;
t_me=mean(t_old);

%variables for project 2
t_new=zeros(377,1);
error_limit=0.01;
time_error_limit=0.01;
rho_chip=7800;
rho_box=2000;
cp_chip=465;
cp_box=840;
t_initial=20;
alpha=k_c/(rho_chip*cp_chip);
dt=60*(dx^2/(4*alpha));
t_error=1;
time_real=0;
time_cal=0;
dt_max=7;
flag_2=0;

while t_error>time_error_limit
    error=1;
    flag_1=0;
    count=0;
    t_n=zeros(377,1)+t_initial;
    while error>error_limit
        %% the code
        for row=1:13
            for col=0:28
                n=row+(col*13);
                if (row==1) %bottom row nodes
                   if (col==0) %node 1
                       a(n,n+1)=k_s;
                       a(n,n+13)=k_s;
                       a(n,n)=-(2*k_s+rho_box*cp_box*dx*dy/(2*dt));
                       const(n,1)=-t_n(n,1)*rho_box*cp_box*dx*dy/(2*dt);
                   end
                   if (col==28) %node 365
                       a(n,n+1)=k_s;
                       a(n,n-13)=k_s;
                       a(n,n)=-(2*k_s+rho_box*cp_box*dx*dy/(2*dt));
                       const(n,1)=-t_n(n,1)*rho_box*cp_box*dx*dy/(2*dt);
                   end
                   if (col>0) && (col<28) %nodes 14:13:352
                       if (n==53) || (n==313) %nodes 53&&313 (connected to the base)
                           a(n,n)=1;
                           const(n,1)=t_base;
                       else 
                           a(n,n+1)=2*k_s;
                           a(n,n+13)=k_s;
                           a(n,n-13)=k_s;
                           a(n,n)=-(4*k_s+rho_box*cp_box*dx*dy/dt);
                           const(n,1)=-t_n(n,1)*rho_box*cp_box*dx*dy/dt;
                       end
                   end
                end
                if ((col==0) || (col==28)) && ((row>1)&&(row<13))
                    if (col==0) %node 2:12
                        a(n,n+1)=k_s;
                        a(n,n-1)=k_s;
                        a(n,n+13)=2*k_s;
                        a(n,n)=-(4*k_s+rho_box*cp_box*dx*dy/dt);
                        const(n,1)=-t_n(n,1)*rho_box*cp_box*dx*dy/dt;
                    else %node 366:376
                        a(n,n+1)=k_s;
                        a(n,n-1)=k_s;
                        a(n,n-13)=2*k_s;
                        a(n,n)=-(4*k_s+rho_box*cp_box*dx*dy/dt);
                        const(n,1)=-t_n(n,1)*rho_box*cp_box*dx*dy/dt;
                    end
                end
                if (row>1) && (row<13)
                   if (row>=9) && ((col>8) && (col<20))
                      if col==9 || col==19
                          if row==9
                              if n==126 %node 126
                                a(n,n+1)=k_average;
                                a(n,n-1)=k_s;
                                a(n,n+13)=k_average;
                                a(n,n-13)=k_s;
                                a(n,n)=-(2*k_average+2*k_s+3*rho_box*cp_box*dx*dy/(4*dt)+rho_chip*cp_chip*dx*dy/(4*dt));
                                const(n,1)=-(qdot*dx*dy/4+t_n(n,1)*(3*rho_box*cp_box*dx*dy/(4*dt)+rho_chip*cp_chip*dx*dy/(4*dt)));
                              else %node 256
                                a(n,n+1)=k_average;
                                a(n,n-1)=k_s;
                                a(n,n-13)=k_average;
                                a(n,n+13)=k_s;
                                a(n,n)=-(2*k_average+2*k_s+3*rho_box*cp_box*dx*dy/(4*dt)+rho_chip*cp_chip*dx*dy/(4*dt));
                                const(n,1)=-qdot*dx*dy/4-t_n(n,1)*(3*rho_box*cp_box*dx*dy/(4*dt)+rho_chip*cp_chip*dx*dy/(4*dt));
                              end
                          else
                              if col==9 %node 127:129
                                a(n,n+1)=k_average;
                                a(n,n-1)=k_average;
                                a(n,n+13)=k_c;
                                a(n,n-13)=k_s;
                                a(n,n)=-(k_s+k_c+2*k_average+rho_box*cp_box*dx*dy/(2*dt)+rho_chip*cp_chip*dx*dy/(2*dt));
                                const(n,1)=-(qdot*dx*dy/2)-(t_n(n,1)*(rho_box*cp_box*dx*dy/(2*dt)+rho_chip*cp_chip*dx*dy/(2*dt)));
                              else %257:259
                                a(n,n+1)=k_average;
                                a(n,n-1)=k_average;
                                a(n,n+13)=k_s;
                                a(n,n-13)=k_c;
                                a(n,n)=-(k_c+k_s+2*k_average+rho_box*cp_box*dx*dy/(2*dt)+rho_chip*cp_chip*dx*dy/(2*dt));
                                const(n,1)=-(qdot*dx*dy/2)-(t_n(n,1)*(rho_box*cp_box*dx*dy/(2*dt)+rho_chip*cp_chip*dx*dy/(2*dt)));
                              end
                          end
                      else
                          if (row == 9) %node 139:13:243
                            a(n,n+1)=k_c;
                            a(n,n-1)=k_s;
                            a(n,n+13)=k_average;
                            a(n,n-13)=k_average;
                            a(n,n)=-(k_s+k_c+2*k_average+rho_box*cp_box*dx*dy/(2*dt)+rho_chip*cp_chip*dx*dy/(2*dt));
                            const(n,1)=-(qdot*dx*dy/2)-(t_n(n,1)*(rho_box*cp_box*dx*dy/(2*dt)+rho_chip*cp_chip*dx*dy/(2*dt)));
                          else %on the chip
                            a(n,n+1)=k_c;
                            a(n,n-1)=k_c;
                            a(n,n+13)=k_c;
                            a(n,n-13)=k_c;
                            a(n,n)=-(4*k_c+rho_chip*cp_chip*dx*dy/dt);
                            const(n,1)=-qdot*dx*dy-t_n(n,1)*(rho_chip*cp_chip*dx*dy/dt);
                          end
                      end
                   else %on the box
                       if (col>0) && (col<28)
                          a(n,n+1)=k_s;
                          a(n,n-1)=k_s;
                          a(n,n+13)=k_s;
                          a(n,n-13)=k_s;
                          a(n,n)=-(4*k_s+rho_box*cp_box*dx*dy/dt);
                          const(n,1)=-t_n(n,1)*(rho_box*cp_box*dx*dy/dt);
                       end
                   end
                end
                if (row==13) %top row nodes
                    if (col==0) %node 13
                        a(n,n-1)=k_s;
                        a(n,n+13)=k_s;
                        a(n,n)=-(2*k_s+h_inf*dx+rho_box*cp_box*dx*dy/(2*dt));
                        const(n,1)=-h_inf*dx*t_inf-t_n(n,1)*rho_box*cp_box*dx*dy/(2*dt);
                    end
                    if (col==28) %node 377
                        a(n,n-1)=k_s;
                        a(n,n-13)=k_s;
                        a(n,n)=-(2*k_s+h_inf*dx+rho_box*cp_box*dx*dy/(2*dt));
                        const(n,1)=-h_inf*dx*t_inf-t_n(n,1)*rho_box*cp_box*dx*dy/(2*dt);
                    end
                    if (col>0)&&(col<28)
                        if col==9  %node 130
                            a(n,n-1)=2*k_average;
                            a(n,n+13)=k_c;
                            a(n,n-13)=k_s;
                            a(n,n)=-(2*h_inf*dx+k_s+k_c+2*k_average+(rho_box*cp_box+rho_chip*cp_chip)*dx*dy/(2*dt));
                            const(n,1)=-(2*h_inf*dx*t_inf+qdot*dx*dy/2+t_n(n,1)*(rho_box*cp_box+rho_chip*cp_chip)*dx*dy/(2*dt));
                        end
                        if col==19 %node 260
                            a(n,n-1)=2*k_average;
                            a(n,n+13)=k_s;
                            a(n,n-13)=k_c;
                            a(n,n)=-(2*h_inf*dx+k_s+k_c+2*k_average+(rho_box*cp_box+rho_chip*cp_chip)*dx*dy/(2*dt));
                            const(n,1)=-(2*h_inf*dx*t_inf+qdot*dx*dy/2+t_n(n,1)*(rho_box*cp_box+rho_chip*cp_chip)*dx*dy/(2*dt));
                        end
                        if (col>9) && (col<19) %node 143:13:247
                            a(n,n-1)=2*k_c;
                            a(n,n+13)=k_c;
                            a(n,n-13)=k_c;
                            a(n,n)=-(4*k_c+2*h_inf*dx+rho_chip*cp_chip*dx*dy/dt);
                            const(n,1)=-2*h_inf*dx*t_inf-qdot*dx*dy-t_n(n,1)*(rho_chip*cp_chip*dx*dy/dt);
                        end
                        if (col<9)||(col>19) %node 26:13:117 && node 273:13:364
                            a(n,n-1)=2*k_s;
                            a(n,n+13)=k_s;
                            a(n,n-13)=k_s;
                            a(n,n)=-(4*k_s+2*h_inf*dx+rho_box*cp_box*dx*dy/dt);
                            const(n,1)=-2*h_inf*dx*t_inf-t_n(n,1)*(rho_box*cp_box*dx*dy/dt);
                        end
                    end
                end
            end
        end
        %% finishing works
        count=count+1;
        t_new=a\const;
        t_average=mean(t_new);
        error=(t_me-t_average)/t_me;
        t_n=t_new;
        if t_n(193,1)>= t_limit && flag_1==0
           time_finish=count*dt;
           flag_1=1;
        end
    end
    time_real=count*dt;
    t_error=(time_real-time_cal)/time_real;
    if flag_2==0
        time_cal=time_real;
         flag_2=1;
         dt=dt_max;
     end
     if flag_2==2
         dt=dt-0.1;
     end
     flag_2=2;
end

dx=1:29;
dy=1:13;
[xt,yt]=meshgrid(dx,dy);
z=t_n(yt+((xt-1)*13));
surfc(xt,yt,z);
xlabel('Column - [1-29]');
ylabel('Row - [1-13]');
zlabel('Tempreture - [C]');
colorbar;
display(dt_max);
X=['Calculated time is: ', num2str(floor(time_cal/60)),' minutes and ', num2str(time_cal-(floor(time_cal/60)*60)),' seconds'];
disp(X);
X=['Real time is: ', num2str(floor(time_real/60)),' minutes and ', num2str(time_real-(floor(time_real/60)*60)),' seconds'];
disp(X);
X=['Elapsed time until the chip reaches the temperature of ', num2str(t_limit) ,' Celsius is: ', num2str(floor(time_finish/60)),' minutes and ', num2str(time_finish-(floor(time_finish/60)*60)),' seconds'];
disp(X);
