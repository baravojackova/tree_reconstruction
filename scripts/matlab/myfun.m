%% my functions

classdef myfun 
    methods (Static = true)
        function  info_res = get_result (QSMs, var)
            % create table
            n = length (QSMs);
            n_v = length (var);
            % create variable types
            for i=1:n_v
                var_typ (i) =  {'double'};
            end
            var_typ(1) = {'string'};
            %
            info_res = table ('Size',[n n_v],'VariableTypes',var_typ, ...
            'VariableNames', var);
            for j=1:n
               for i=1:n_v
                 columnname = var {i};
                 if strcmp(columnname, 'name') || strcmp(columnname, 'tree')
                    info_res.(columnname)(j) = QSMs(j).rundata.inputs.(columnname);  
                    else
                    info_res.(columnname)(j) = QSMs(j).treedata.(columnname);  
                 end
               end
            end
        end
        %
        function geom_all_step = result_ansys (QSMs,n_opt)
            geom_all_step = cell (n_opt,1);
            for i=1:n_opt
                geom_all_step{i,1}=zeros(length(QSMs(i).cylinder.length)+1,11);
                geom_all_step{i,1}(1,:)=[0 1 2 3 4 5 6 7 8 9 10];
                geom_all_step{i,1}(1:end,1)=(0:length(QSMs(i).cylinder.length));
                geom_all_step{i,1}(2:end,2:4)=QSMs(i).cylinder.start;
                geom_all_step{i,1}(2:end,5)=QSMs(i).cylinder.length;
                geom_all_step{i,1}(2:end,6)=QSMs(i).cylinder.radius;
                geom_all_step{i,1}(2:end,7:9)=QSMs(i).cylinder.axis;
                geom_all_step{i,1}(2:end,10)=QSMs(i).cylinder.parent;
                geom_all_step{i,1}(2:end,11)=QSMs(i).cylinder.BranchOrder;
            end
        end
    end
end
