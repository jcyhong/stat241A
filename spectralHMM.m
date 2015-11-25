data = Quandl.get('YAHOO/INDEX_GSPC','start_date','2010-11-23','end_date','2015-11-23','transformation','rdiff');
Close = data.Close;
obs = Close.data;

% pulls in SP500 data in using quandl

tenday=zeros(125,1);
for i=1:125
    j=10*i-9;
    tenday(i)=prod(1+obs(j:(j+9)));
end

tenday = tenday - 1;
% the above code generates a list of rates of increase for 10 day
% intervals. do this to avoid having really small observation values

timestep = (1:length(tenday))';
hmm = [timestep HMM_obs(tenday,6)];
tenday = [timestep tenday];

%create empirical estimates of P1, P21, P3_x_1
