%% Assign Weight
%weight_var = -1/mean(AverageVariability);
%weight_demand = 1/mean(Avg_Traffic);
weight_var_multiplier = 1
weight_demand_multiplier = 1


%% Set Parameters
max_hours_per_week = 8; % Maximum number of hours a driver can work per week
time_slot_available = ones(7,24); % Represent the timeslot when the driver is available
%time_slot_available(2,:) = 0;
%time_slot_available(4,:) = 0;
%time_slot_available(6:7,:) = 0;
%time_slot_available(:,1:18) = 0;

region_avilable = ones(5,1); % Indicate which of the 5 regions the driver is avilable to go
% Store datasets
%Weekday_Copy = Weekday;
%Timeslot_Copy = Timeslot;
%Avg_Var_Copy = AverageVariability;
%Avg_Demand_Copy = Avg_Traffic;
%Region_Copy = Region;

% Restore datasets
Weekday = Weekday_Copy;
Timeslot = Timeslot_Copy;
AverageVariability = Avg_Var_Copy;
Avg_Traffic = Avg_Demand_Copy;
Region = Region_Copy;


%% Import Data
% Right click on the "Data_Demand_Variability.csv", select "Import Data"

%% Adjust the data based on the given parameters - driver's preference
weekday_choices = 1:7;
time_choices = [0 1 7 8 16 17 18 19 20 21 22 23];
timeslot_choices = [0 1 2 3 4 5]; % Timeslots that cooresponds to the above time_choices
n_weekdays_choices = size(weekday_choices,2);
n_timeslot_choices = size(weekday_choices,2);

% Initialize indices of element to be removed: all are 0 in the beginning
indices_to_be_rm = zeros(size(Weekday));
for d = weekday_choices
    for t = time_choices
        if time_slot_available(d,t+1) == 0 % Mark as in-available
           % Find the timeslot and weekday that cooresponding to that hour
           t_ind = find(time_choices==t);
           if mod(t_ind,2) == 0
               timeslot_to_be_rm = timeslot_choices(floor(find(time_choices==t)/2));
           else
               timeslot_to_be_rm = timeslot_choices(floor(find(time_choices==t)/2)+1);
           end
           weekday_to_be_rm = d;
           % Mark the indices of the elements to be removed
           indices_to_be_rm(Weekday == weekday_to_be_rm & Timeslot == timeslot_to_be_rm) = 1;
        end
    end
end

% Mark unavailable reigons
% We need to add 1 to the index of region
if sum(region_avilable)<5
    indices_to_be_rm(Region == find(region_avilable==0)-1) = 1;
end

% Removed all marked instances
% We need to backwards
if sum(indices_to_be_rm) >0
    ind_to_be_removed = find(indices_to_be_rm == 1);
    for i = sum(indices_to_be_rm):-1:1
        ind = ind_to_be_removed(i);
        Weekday(ind) = [];
        Timeslot(ind) = [];
        AverageVariability(ind) = [];
        Avg_Traffic(ind) = [];
        Region(ind) = [];
    end
end

% Number of variables
n_x = size(AverageVariability,1);
%% Formulate the IP
weight_var = -1/mean(AverageVariability)*weight_var_multiplier;
weight_demand = 1/mean(Avg_Traffic)*weight_demand_multiplier;
%  A vector of cost coefficients
f = transpose(weight_var.*AverageVariability + weight_demand.*Avg_Traffic);
% Integer Variables
intcon = 1:n_x;

% Constraint 1: Maximum number of hours per week
A1 = ones(1,n_x); % Lhs
b1 = max_hours_per_week/2;  % rhs

% Constraint 2: You cannot be in 2 regions during the same time
% Initialize A2 and b2
A2 = zeros(n_weekdays_choices*n_timeslot_choices,n_x);
b2 = ones(n_weekdays_choices*n_timeslot_choices,1); % Upper limit should be 1

row_count = 0;
for d = weekday_choices
    for t = timeslot_choices
        row_count = row_count+1;
        ind_marked = (Weekday == d & Timeslot == t);
        A2(row_count,ind_marked) =1;
    end
end
    
A = [A1;A2];
b = [b1;b2];
lb = zeros(n_x,1); 
ub = ones(n_x,1);
[x,fval] = intlinprog(-f,intcon,A,b,[],[],lb,ub);

%% Present the Solution
select_list = find(x==1);

% Set a basis for the weekday
day_number_benchmark = datetime(2017,04,9);

% Set a list of our pre-defined timeslot
timeslots = {'12AM - 2AM'; '7AM - 9AM'; '4PM - 6PM'; '6PM-8PM'; '8PM-10PM'; '10PM-12AM'};

for i = 1:sum(x)
    select_time = Timeslot(select_list(i));
    select_day = Weekday(select_list(i));
    select_region = Region(select_list(i));
    
    weekday_name = day(day_number_benchmark+select_day,'name');
    timeslot = timeslots(select_time +1);
    string_output = ['Recommendation ' num2str(i) ': Region ' num2str(select_region) ', ' timeslot ', ' weekday_name];
    disp(string_output)
end
%% Example IP Solver - Minimize the cost
%  vector of cost coefficients
%f = [350*5,330*3,310*4,280*6,500,450,400,100];
% Integer Variables
%intcon = 1:4;
%Aeq = [5,3,4,6,1,1,1,1;
%    5*0.05,3*0.04,4*0.05,6*0.03,0.08,0.07,0.06,0.03;
%    5*0.03,3*0.03,4*0.04,6*0.04,0.06,0.07,0.08,0.09];
%beq = [25;1.25;1.25];
%lb = zeros(8,1);
%ub = ones(8,1);
%ub(5:end) = Inf; % No upper bound on noninteger variables
%[x,fval] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
%[x,fval] = intlinprog(f,intcon,[],[],Aeq,beq,lb,ub);
