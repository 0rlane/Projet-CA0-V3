load ps.pts
load Omega.pts
load NM.pts

N=ps(1,1);
ps2=sortrows(ps([2:end],:));
tri=delaunay(ps2(:,2),ps2(:,3));
hold on;
trimesh(tri,ps2(:,2),ps2(:,3),zeros(size(ps2(:,2))));

% Affichage des centres inscrit des triangles
plot(Omega(:,1),Omega(:,2),'ro')

% Affichage des Mi
plot([NM(:,1);NM(:,3);NM(:,5)],[NM(:,2);NM(:,4);NM(:,6)],'b+')

hold off;

% Generation d'un fichier de doubles contenant les corrdonnees des points
% dans l'odre des sommets
ficpoints='points.pts';
[fid, message]=fopen(ficpoints,'w');
if (fid < 0) error([message,'(fichier ',ficpoints,')']),end
fprintf(fid,'%d\n',N);
fprintf(fid,'%f %f %f\n',ps2');
fclose(fid);


% Generation d'un fichier d'entiers contenant la triangulation
fictri = 'listri.dat';
[fid, message] = fopen(fictri,'w');
if (fid < 0) error([message,' (fichier ',fictri,')']), end
fprintf(fid,'%d\n',length(tri));
fprintf(fid,'%d %d %d\n',tri');
fclose(fid);
