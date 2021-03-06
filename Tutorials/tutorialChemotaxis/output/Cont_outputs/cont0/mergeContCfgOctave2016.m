#! /usr/bin/octave -qf

clear
filepath = './'
info = load ([filepath 'cfg-lings/c.info'])
nRank = info(1);
nMol = info(2);
cfg_intv = info(3);
stopAt  = info(4);
entry_count = info(5);
startAt = 0;

%~ tick = 100;
%%%% indicate whether you want to remove cfg-lings as you build Cfg files
remove_files = 0;
for tick = startAt:cfg_intv:stopAt
	try 
		table0 = [];
		c_type_file0 = [filepath 'cfg-lings/c.0_type'] ;
		for rank = 0:nRank-1
			cfglings_file = [filepath 'cfg-lings/c.' int2str(tick) '_' int2str(rank) ];
			try
				cfg0 = load (cfglings_file);
			catch
				cfg0 = [];
			end
			%~ [m , n] = size(cfg0);
			
			%~ for i = 1 : m
				%~ table0(cfg0(i,1)+1,:) = cfg0(i,2:n);
			%~ end
			table0 = [table0 ; cfg0];
		end
		c_head_file = [filepath 'cfg-lings/c.0_head'] ;% Edited By Vi on 2016/04/07
		%~ c_head_file = [filepath 'cfg-lings/c.' int2str(tick) '_head'] ;
		c_body_file = [filepath 'cfg-lings/c.' int2str(tick) '_body'] ;
		cfg_final_file = [filepath 'Cfg/s.' int2str(tick)]

		dlmwrite(c_body_file, table0, ' ');
		linux_command0 = [ 'rm -f ' cfg_final_file  ];
		
		system(linux_command0);
		
		fileID = fopen(cfg_final_file,'w');
		fprintf(fileID,'Number of particles = %d\n', nMol);
				
		fclose(fileID);
		
		linux_command1 = [ 'cat ' c_head_file ' >> ' cfg_final_file];
		linux_command1b = [ 'cat ' c_type_file0 ' >> ' cfg_final_file];
		linux_command2 = [ 'cat ' c_body_file ' >> ' cfg_final_file];


		system(linux_command1);
		system(linux_command1b);
		system(linux_command2);

		linux_command_rmbody = [  'rm -f ' c_body_file ];

		system(linux_command_rmbody);


		if remove_files
			cfglings_file_list = [];
			for rank = 0:nRank-1
				cfglings_file = ['../../cfg-lings/c.' int2str(tick) '_' int2str(rank) ];
				cfglings_file_list = [cfglings_file_list ' ' cfglings_file ] ;
			end 
			linux_command_rm = [ 'rm -f ' cfglings_file_list ' ' c_head_file ' ' c_body_file];
			system(linux_command_rm);

		end
	catch
	%~ 'File is currently in use by another process or thread'
	end
end
