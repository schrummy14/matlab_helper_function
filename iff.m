function y = iff(condition,t_condition,f_condition)

    y = condition.*t_condition + (1-condition).*f_condition;

end