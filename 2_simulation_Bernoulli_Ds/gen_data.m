%% D
Draw=readtable('palmer_scratch/final_analysis_Dfs/Stratified_nat_13_D.csv');
Draw=table2array(Draw(:,2:9019));

p=readtable('Unified/final_analysis/output/Stratified_nat_13_p.csv');
p=table2array(p(:,2));
D13_nat = diag(1./p) * Draw * diag(1./p);
[eigen_vec,eigen_value]=eigs(D13_nat,1);

%%
D=D13_nat;

%% 
N=4509;
X=readtable('Unified/final_analysis/output/x.csv');
X=table2array(X(:,2:8));
intercept1=repmat([1,0],1,N);
intercept3=repmat([0,1],1,N);
X=repelem(X,2,1);
X=[intercept1',intercept3',X];

%% 
filter= (-1).^(1:(2*N));
filter=filter';
X=diag(filter) * X;

%%
XDX = X' * D * X/N;

Optimal_Project = eye(2*N) - X * inv(XDX) * X' * D/N;

Optimal_Sandwich = Optimal_Project' * D * Optimal_Project;

%% 
XX= (diag(filter) * X)' * ((diag(filter) * X))/N;

WLS_Project= eye(2*N) - X * inv(XX) * X'/N;

WLS_Sandwich = WLS_Project' * D * WLS_Project;


%% 
[eigen_vec,eigen_value] = eigs(Optimal_Sandwich/N-WLS_Sandwich/N,20);

%% 
[eigen_vec,eigen_value] = eigs(WLS_Sandwich/N-D/N,20);

%%
for i=9:20
    i
    eigen_value(i,i)
    y=eigen_vec(:,i);
    y = round(0.5+diag(filter)*y);

    HT_var_new = (y.*filter)' *  D * (y.*filter) /N;
    WLS_var_new = (y.*filter)' *  WLS_Sandwich * (y.*filter) /N;
    Optimal_var_new = (y.*filter)' *  Optimal_Sandwich * (y.*filter) /N;
    (WLS_var_new-Optimal_var_new)/ Optimal_var_new
    (HT_var_new-WLS_var_new)/ WLS_var_new

end


%%
 pred=X * inv(XX) * X'/N * (y.*filter);

 alpha= inv(pred' * D * pred/N)* pred' * D*(y.*filter)/N;

 %%
WLS_Project_a= eye(2*N) - alpha* X * inv(XX) * X'/N;

WLS_Sandwich_a = WLS_Project_a' * D * WLS_Project_a;

WLS_var_new_a = (y.*filter)' *  WLS_Sandwich_a * (y.*filter) /N;


%%
HT_var = (y.*filter)' *  D * (y.*filter) /N;
%%

y=sum(eigen_vec(:,1:i),2);
y=round(0.5+diag(filter)*y);
csvwrite('/vast/palmer/home.grace/hc654/Unified/final_analysis/13.csv',y);