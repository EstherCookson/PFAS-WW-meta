function output_cell = cell_str_2_num(input_cell)

    for i = 1:length(input_cell)
        
        str_val = input_cell{i};
        
        if contains(str_val,'<')
            str_val2 = str2double(str_val(2:end))* 0.5; % *rand(1,1)
            input_cell{i}=str_val2;

        else
            input_cell{i}=str2double(str_val);
        end 
    end 

	output_cell = cell2mat(input_cell);

end