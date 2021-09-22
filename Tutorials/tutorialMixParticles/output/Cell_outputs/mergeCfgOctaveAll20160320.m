#! /usr/bin/octave -qf

clear
%~ data = ['./cell0/';'./cell1/';'./cell2/'];
%~ filepath = cellstr(data)
%~ filepath = celldata{0+1}
%~ filepath1 = './cell1/';
out_filepath = './'
info{1+0} = load ('./cell0/cfg-lings/c.info')



nRank = info{1+0}(1)
cfg_intv = info{1+0}(3)
stopAt  = info{1+0}(4)
%~ entry_count = info{1+0}(5)
total_number_of_cells = info{1+0}(5)
startAt = 0;

for cell_count_n = 0 : total_number_of_cells - 1
	filepath{1+cell_count_n} = ['./cell' int2str(cell_count_n) '/']
	info{1+cell_count_n} = load ([filepath{1+cell_count_n} 'cfg-lings/c.info'])
end
%~ nMol{1+0} = info{1+0}(2)
%~ nMol{1+1} = info{1+1}(2)

%~ tick = 100;
%%%% indicate whether you want to remove cfg-lings as you build Cfg files
remove_files = 0;
for tick = startAt:cfg_intv:stopAt
	try 
		c_head_file = [filepath{1+0} 'cfg-lings/c.' int2str(tick) '_head'] ;
		cfg_final_file = [out_filepath 'Cfg/s.' int2str(tick)]
		
		linux_command_cfg_final_file = [ 'rm -f ' cfg_final_file  ];
		system(linux_command_cfg_final_file);
		
		
		
		
		nMolTotal = 0;
		for cell_count_n = 0 : total_number_of_cells - 1
			nMolTotal = nMolTotal + info{1+cell_count_n}(2);
		end
		fileID = fopen(cfg_final_file,'w');
		fprintf(fileID,'Number of particles = %d\n', nMolTotal);
		fclose(fileID);
		
		
		
		linux_command1 = [ 'cat ' c_head_file ' >> ' cfg_final_file];
		system(linux_command1);
		%% 0000000000000000000000000000000000000000000000000000000000000
		
		for cell_count_n = 0 : total_number_of_cells - 1
			if info{1+cell_count_n}(2) > 0
				table0 =[];
				c_type_file{1+cell_count_n} = [filepath{1+cell_count_n} 'cfg-lings/c.0_type'] ;
				for rank = 0:nRank-1
					cfglings_file = [filepath{1+cell_count_n} 'cfg-lings/c.' int2str(tick) '_' int2str(rank) ];
					try
						cfg0 = load (cfglings_file);
					catch
						cfg0 = [];
					end
					[m , n] = size(cfg0);
					
					for i = 1 : m
						table0(cfg0(i,1)+1,:) = cfg0(i,2:n);
					end
				end
				c_body_file{1+cell_count_n} = [filepath{1+cell_count_n} 'cfg-lings/c.' int2str(tick) '_body'] ;
				%~ cfg_final_file = [filepath 'Cfg/s.' int2str(tick)]

				dlmwrite(c_body_file{1+cell_count_n}, table0, ' ');
				
				
				
				
				linux_command1b = [ 'cat ' c_type_file{1+cell_count_n} ' >> ' cfg_final_file];
				linux_command2 = [ 'cat ' c_body_file{1+cell_count_n} ' >> ' cfg_final_file];
				linux_command_rmbody = [  'rm -f ' c_body_file{1+cell_count_n} ];
				system(linux_command1b);
				system(linux_command2);
				system(linux_command_rmbody);
			end
		end

		%~ cell_count_n = 1;
		%~ 
%~ 
		%~ 
		%~ linux_command1b = [ 'cat ' c_type_file{1+cell_count_n} ' >> ' cfg_final_file];
		%~ linux_command2 = [ 'cat ' c_body_file{1+cell_count_n} ' >> ' cfg_final_file];
		%~ linux_command_rmbody = [  'rm -f ' c_body_file{1+cell_count_n} ];
		%~ system(linux_command1b);
		%~ system(linux_command2);
		%~ system(linux_command_rmbody);
		
		

		if remove_files
			cfglings_file_list = [];
			for rank = 0:nRank-1
				cfglings_file = ['../../cfg-lings/c.' int2str(tick) '_' int2str(rank) ];
				cfglings_file_list = [cfglings_file_list ' ' cfglings_file ] ;
			end 
			linux_command_rm = [ 'rm -f ' cfglings_file_list ' ' c_head_file ' ' c_body_file{1+0}];
			system(linux_command_rm);

		end
	catch
	%~ 'File is currently in use by another process or thread'
	end
end
