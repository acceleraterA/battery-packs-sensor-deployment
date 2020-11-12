function[x,y]=num2node(num,m,n)
y=mod(num,n);%index of column equals to resudual after devided by total number of rows 
x=ceil(num/n);% index of row equals to num/number of column 
end
