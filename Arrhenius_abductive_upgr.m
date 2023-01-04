% The code is designed to find the expression known as Arrhenius equation. The operators, constants and blacklist are chosen
% accordingly.


%%% User provided inputs and hyperparameters%%% 

%mathematical_constants = {'pi','exp(1)','double(eulergamma)','((1+sqrt(5))/2)','1','2','3','4','-1'}; %Pi, Euler's number, Euler's constant, and Golden ratio
mathematical_constants = {'exp(1)'};

%fundamental_physical_constants = {'physconst(''LightSpeed'')', 'physconst(''Boltzmann'')','6.63e-34','2.898e-3','6.6738e-11'}; %Speed_of_light, Boltzmann constant, Planck_constant, Wein's Law constant, Universal gravitational constant
fundamental_physical_constants = {'8.31446261815324'};

variables = {'A'};

%operators = {'+','-','./','.*','.^','sqrt(','(',')','log(','ln(','sin(','cos(','tan(','ctg(', 'factorial('};
operators = {'.^',')','(','/','.*'};

nonfundamental_physical_constants = {'x(1)'}; % preserve this sintax because of the optimization input format requirements

blacklist = {'exp(1)e','exp(1)p','exp(1)2','exp(1)A','exp(1)(','piexp(1)','pipi','pi2','piA','pi(','2exp(1)','2pi','22','2A','2(','Aexp(1)',...
    'Api','A2','A(','.^.^','.^.*','.^)','.^./','.*.^','.*.*','.*)','.*./','(.^','(.*','()','(./',')exp(1)',')pi',')2',')A',')(','./.^','./.*','./)','././','x(1)x(2)','x(2)x(1)'};
    

% for factorial calculation we chose three natural numbers to be used as
% checkpoints

var_at_checkpoints = [-0.04,-0.0384,-0.0340];

Y_c = [4.323552476185454e+07,4.330475697247247e+07,4.349571770868780e+07];%.*[unifrnd(0.99,1.01),unifrnd(0.99,1.01),unifrnd(0.99,1.01)];

k = 5; %the length of the candidate_law_syntax

tolerance = 0.0001;




%%% Main code %%%


% create struct datatype for optimization purposes
for i=1:size(variables,2)
    eval(['var.',variables{1,i},'=',mat2str(var_at_checkpoints(i,:)),';']);
end

m = size(Y_c,2);

% Concatenate three lists into one
elements = [mathematical_constants, fundamental_physical_constants, variables, operators, nonfundamental_physical_constants];

% Create variations with repetitions using all elements
n = size(elements,2); % total number of basic strings that can be used to construct the law syntax
no_of_patterns = n^k;


c=1;
solutions={};
nonfundamental_constants = {};
best_candidate = {};
distance_to_truth = inf;
generator = nextstring(n,k);
progressbar
x=nan;
for i=1:no_of_patterns
    indices = generator();
    string_subset = elements(1,indices);
    candidate = strjoin(string_subset,'');
    
    if ~contains(candidate, blacklist)
        
        no_of_nonfundamental_constants = sum(ismember(unique(string_subset),nonfundamental_physical_constants));
        
        if no_of_nonfundamental_constants>0
            
            % specifiy starting position in optimization
            x0 = ones(1,no_of_nonfundamental_constants);
  
            try
                % replace variable symbols with var.variable strings in
                % candidate expression
                upgr_candidate = candidate;
                for k=1:size(variables,2)
                    eval(['upgr_candidate = replace(upgr_candidate,''',variables{1,k},''',''var.',variables{1,k},''');']);
                end

                % define anonymous function
                eval(['non_fund_fitting = @(x,xdata) ', upgr_candidate, ';']);
                
                % fit the function into candidate law
                opts = optimset('Display','off');
                eval(['x = fminunc(@(x) sum((Y_c - ',upgr_candidate,').^2),x0,opts);'])
            catch
                progressbar(i/no_of_patterns)
                continue
            end
        end

        % define flags to track the success of validations at every checkpoint
        flags=zeros(1,m);
        value_accu = zeros(1,m);
        try 
            for j=1:m
                for k=1:size(variables,2)
                    eval([variables{1,k},'=',mat2str(var_at_checkpoints(k,j)),';']);
                end
                Y=Y_c(1,j);
                current_value = eval(candidate);
                value_accu(1,j)=current_value;

                if abs((current_value-Y)/Y) < tolerance
                    flags(1,j) = 1;
                end
            end  
            if sum(flags)==m
                solutions{c,1}=candidate;
                nonfundamental_constants{c,1}=x;
                if sum(abs(value_accu-Y_c))<distance_to_truth
                    distance_to_truth = sum(abs(value_accu-Y_c));
                    best_candidate = candidate;
                end
                x=nan;
                c=c+1;
            end
        catch
            progressbar(i/no_of_patterns)
            continue
        end
    end
    progressbar(i/no_of_patterns)
end
    
    
    


