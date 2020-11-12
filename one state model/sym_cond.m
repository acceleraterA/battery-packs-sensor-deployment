m = input('Please enter the number of the row of the battery system:'); 
n = input('Please enter the number of the column of the battery system:');
if n==1%1*n and m*1 are the same for string of batteries
    n=m;
    m=1;
end

A = sym(zeros(m,n));
syms Rc Ru Rc Cs Rcc Cf Cc

for i=2:2:2*m*n
    %%Ts
    %cell to cell conduction
    if m>1
        A(i,i)=-(2/(Ru*Cs)+4/(Rcc*Cs));%all cells except the cells on the boundary have 4 Rccs -4.6546
        A(2,2)=-(2/(Ru*Cs)+2/(Rcc*Cs));
        A(2*n,2*n)=-(2/(Ru*Cs)+2/(Rcc*Cs));
        A(2*(m-1)*n+2,2*(m-1)*n+2)=-(2/(Ru*Cs)+2/(Rcc*Cs));
        A(2*m*n,2*m*n)=-(2/(Ru*Cs)+2/(Rcc*Cs));%4 cells on the corner with 2 Rccs    -2.4323
        if n>2
            for j =4:2:2*(n-1)
                A(j,j)=-(2/(Ru*Cs)+3/(Rcc*Cs));%top row except the corner -3.5435
            end
            for j=2*(m-1)*n+4:2:2*m*n-2
                A(j,j)=-(2/(Ru*Cs)+3/(Rcc*Cs));%bottom row except the corner
            end
        end
        if m>2
            for j=2*n+2:2*n:2*((m-2)*n+1)
                A(j,j)=-(2/(Ru*Cs)+3/(Rcc*Cs));%first column except the corner
            end
            for j=(2*2*n):2*n:(m-1)*2*n
                A(j,j)=-(2/(Ru*Cs)+3/(Rcc*Cs));%last column except corner, cells on the boundary with 3 Rccs
            end
        end
        if mod(i/2,n)~=0
            A(i,i+2)=1/(Rcc*Cs);%right connect rcc except the last column 1.1111
        end
        if mod(i/2,n)~=1
            A(i,i-2)=1/(Rcc*Cs);%left connect rcc except the first column
        end
        if i>2*n+1
            A(i,i-2*n)=1/(Rcc*Cs);%above connect rcc except the first row
        end
        if i<=2*(m-1)*n
            A(i,i+2*n)=1/(Rcc*Cs);%below connect rcc except the last row
        end

    end
end   
A(1:2:2*m*n,:)=[];
A(:,1:2:2*m*n)=[];
