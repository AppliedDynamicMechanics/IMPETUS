#! /usr/bin/octave -qf
clf
close all

MSD = 'MSD/files/MSD.out';
RDF = 'RDF/files/gr_current.out';

w = 560;
l = 840;
%~ l = 420;
figure(1)
set (figure(1),'Position', [300 100 2*w 2*l] )

while(true)

try 
	a = load(MSD); 
	subplot(2,1,1) ;
	%~ subplot(2,1,1,'FontSize',14,'FontWeight','bold') ;
	plot(a(:,1),a(:,2), 'linewidth',2);
	grid on;

	xlabel('time','FontSize',16,'FontWeight','bold');
	ylabel('MSD','FontSize',16,'FontWeight','bold');
catch
end


try 
	a = load(RDF); 
	subplot(2,1,2) ;
	%~ subplot(2,1,2,'FontSize',14,'FontWeight','bold') ;
	plot(a(:,1),a(:,2), 'linewidth',2);
	grid on;
	xlabel('radius','FontSize',16,'FontWeight','bold');
	ylabel('G(r)','FontSize',16,'FontWeight','bold');
catch
end



pause(5);


end
