%Reference : Aaron Danner tutorial & Mathematica code (NUS/UIUC), & Kevin D. Welsher FORTRAN
%code
%Author: Muhanad Zaki
%Ref.(Danner): http://www.ece.nus.edu.sg/stfpage/eleadj/pseudopotential.htm
%Ref.(Welsher): http://large.stanford.edu/courses/ap272/welsher1/
%Empirical Pseudopotential Form factors are from Marvin L. Cohen and T.K. Bergstresser,
%?Band Structures and Pseudopotential Form Factors forFourteen Semiconductors of the Diamond and Zinc-blende Structures,?
% Phys. Rev. 141, 2, p. 556, 1966.

clear;
tic
T=[1/8 1/8 1/8];   %Structure factor
%Start generating Brillouin zone vectors
n=(5^3)-1 ;
i=1;
%steps_relative = [2 1 3 2];
steps_relative = [2 1 2 2];
resolution = 10;
steps = steps_relative*resolution;

k_start = [0 0 0]; %Gamma
%k_start = [1/2 1/2 1/2]; 
endpoints = zeros(4,3);

%GaAs
%{
endpoints(1,:) = [0 0 0];
endpoints(2,:) = [1 0 0];
endpoints(3,:) = [3/4 3/4 0];
endpoints(4,:) = [0 0 0];
gamma_endpoint = 1;
%}


endpoints(1,:) = [3/4 3/4 0];
endpoints(2,:) = [1 0 0];
endpoints(3,:) = [0 0 0];
endpoints(4,:) = [1/2 1/2 1/2];
gamma_endpoint = 3;

k_vectors_proper=zeros(3,sum(steps)+1); %k is a vector moving through the Brillouin zone. 

k_counter = 2;
k_vectors_proper(:,1) = k_start(:); %Add starting point
k_current = k_start(:);
gamma = [0 0 0];
if gamma_endpoint == 0
    gamma = [0 0 0];
end;
for i = 1:4
    
    k_end = endpoints(i,:); %Save end point for current path
    k_step = [(k_end(1)-k_current(1))/steps(i) (k_end(2)-k_current(2))/steps(i) (k_end(3)-k_current(3))/steps(i)]; %Step vector for current path
    for m = 0:steps(i)-1
        k_next = [k_current(1)+k_step(1) k_current(2)+k_step(2) k_current(3)+k_step(3)];
        k_vectors_proper(:,k_counter) = k_next;
        k_current = k_next;
        k_counter = k_counter + 1;
    end;
    if i == gamma_endpoint
        gamma = k_counter-1;
    end;
end;
k_vectors = k_vectors_proper;

n = size(k_vectors);
number_of_vectors = n(2);

%List of constants
m0=9.11E-31; %free electron mass in Kg
a=4.3596E-10;   %SiC
%a = 5.6533e-10; %GaAs
q=1.6E-19;
hbarJs=1.054E-34;
hbareV=6.581E-16;

H_matrix_half_size = 62;
H_matrix_size = 2*H_matrix_half_size+1;
energy = zeros(length(k_vectors),H_matrix_size); 
ctr = 0; %For percentage display

for k_iterator=1:number_of_vectors %Iterate through the k-vectors (move through the B-zone
    
    k=k_vectors(:,k_iterator);  
    H_matrix = zeros(H_matrix_size,H_matrix_size);
    
    %Display stuff
    done = k_iterator/41;
    if done>0.9 && ctr<90
       display('90 %');
       ctr = 90;
    elseif done>0.8 && ctr<80
        display('80 %');
        ctr = 80;
    elseif done>0.7 && ctr<70
        display('70 %');
        ctr = 70;
    elseif done>0.6 && ctr<06
        display('60 %');
        ctr = 60;
    elseif done>0.6 && ctr<60
        display('60 %');
        ctr = 60;
    elseif done>0.5 && ctr<50
        display('50 %');
        ctr = 50;
    elseif done>0.4 && ctr<40
        display('40 %');
        ctr = 40;
    elseif done>0.3 && ctr<30
        display('30 %');
        ctr = 30;
    elseif done>0.2 && ctr<20
        display('20 %');
        ctr = 20;
    elseif done>0.1 && ctr<10
        display('10 %');
        ctr = 10;
    end
         
    for i=-H_matrix_half_size:1:H_matrix_half_size %For each point i B-zone (each k), evaluate hamiltonian matrix. 
        for j=-H_matrix_half_size:1:H_matrix_half_size
            
            if i==j ;
                m=i;
                G=Ggen(m);
                DotProduct=dot((G(:)+k(:)),(G(:)+k(:))) ;
                H_matrix(i+H_matrix_half_size+1,j+H_matrix_half_size+1)=((hbarJs*hbareV)/(2*m0))*(abs(DotProduct*(2*pi/a)^2));

            else
                m=i-j;
                G=Ggen(m);
                V=Vpseudo_SiC(G);
                H_matrix(i+H_matrix_half_size+1,j+H_matrix_half_size+1)=V;

            end;

        end;

    end;

    energy(k_iterator,:) = eig(H_matrix);
    %eval(['energy' num2str(k_iterator) '=eig(H_matrix)';]);

end;




fig=figure;
hax=axes;
for energy_level_iterator = 1:8
    
    hold on;
    engergy_adj = energy-energy(gamma,3);
    plot((0:number_of_vectors-1),engergy_adj((1:number_of_vectors),energy_level_iterator));
    
    %SP1=10;
    %SP2=15;
    %SP3=30;
    
    SP1 = steps(1);
    SP2 = SP1 + steps(2);
    SP3 = SP2 + steps(3);
    
    line([SP1 SP1],get(hax,'YLim'),'Color',[0 0 0])
    line([SP2 SP2],get(hax,'YLim'),'Color',[0 0 0])
    line([SP3 SP3],get(hax,'YLim'),'Color',[0 0 0])
end;
hold off;
toc


