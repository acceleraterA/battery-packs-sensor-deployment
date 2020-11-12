function [A,B] = A_conv_make(Ru,Cs,Cf,R_e,m,n)
A = zeros(2*m*n,2*m*n);
B = zeros(2*m*n,2);
for i=2:2:2*m*n

    %cell to cell conduction
    if m>1
        A(i,i)=-2/(Ru*Cs);%all cells except the cells on the boundary have 4 Rccs -4.6546
        A(2,2)=-2/(Ru*Cs);
        A(2*n,2*n)=-2/(Ru*Cs);
        A(2*(m-1)*n+2,2*(m-1)*n+2)=-2/(Ru*Cs);
        A(2*m*n,2*m*n)=-2/(Ru*Cs);%4 cells on the corner with 2 Rccs    -2.4323
        if n>2
            for j =4:2:2*(n-1)
                A(j,j)=-2/(Ru*Cs);%top row except the corner -3.5435
            end
            for j=2*(m-1)*n+4:2:2*m*n-2
                A(j,j)=-2/(Ru*Cs);%bottom row except the corner
            end
        end
        if m>2
            for j=2*n+2:2*n:2*((m-2)*n+1)
                A(j,j)=-2/(Ru*Cs);%first column except the corner
            end
            for j=(2*2*n):2*n:(m-1)*2*n
                A(j,j)=-2/(Ru*Cs);%last column except corner, cells on the boundary with 3 Rccs
            end
        end
        %coolant and battery surface convection
        if i>2*n &&mod(i/2,n)~=1 %start at second row cell and except for first column which tf = tin
            for j =2:2:2*(i/2-(ceil(i/2/n)-1)*n-1)%ceil(i/2/n)=#column
                if i>2*n && i< 2*(m-1)*n +1&& m>2 %except the first and last row
                    A(i,i-j)=A(i,i-j)+2*(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));
                    A(i,i-j-2*n)=A(i,i-j-2*n)+(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));
                    A(i,i-j+2*n)=A(i,i-j+2*n)+(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));
                elseif i> 2*(m-1)*n+2%last row
                    A(i,i-j)=A(i,i-j)+(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs))+(1-1/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));%current row
                    A(i,i-j-2*n)=A(i,i-j-2*n)+(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));%row above
                end
            end
        elseif i<=2*n %first row except first column
            for j=2:2:i-2
                A(i,j)=A(i,j)+(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs))+(1-1/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));%current row
                A(i,j+2*n)=A(i,j+2*n)+(1-2/(Ru*Cf))^(j/2-1)*(1/(Ru^2*Cf*Cs));%row below
            end
        end
        B(i,1)=R_e/Cs;
        B(i,2)=2/(Ru*Cs)*(1-2/(Ru*Cf))^(i/2-(ceil(i/2/n)-1)*n-1);         

    else% for a string of battries
        %conduction between each battries
        A(i,i)=-(1/(Ru*Cs));
        A(2,2)=-(1/(Ru*Cs));
        A(2*n,2*n)=-(1/(Ru*Cs));%first and last cells have one Rcc
       
        %coolant and battery surface convection
        for j = 2:2:i-2
            A(i,j)=A(i,j)+(1-1/(Ru*Cf))^((i-j)/2-1)/(Ru^2*Cf*Cs);
        end
        B(i,1)=R_e/Cs;
        B(i,2)=1/(Ru*Cs)*(1-1/(Ru*Cf))^(i/2-1);

    end
end    
A(1:2:2*m*n,:)=[];
A(:,1:2:2*m*n)=[];
B=B(2:2:2*m*n,1);

