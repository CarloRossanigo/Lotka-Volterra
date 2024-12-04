%%Matlab simply code for Lotka-Volterra model 
%% codice Matlab semplice per visualizzare il modello Lotka-Volterra
%% addpath for file (made by codeblocks) location
addpath '/Users/your_name' %aggiungo il percorso per il file dei risultati

k=load('Lotka'); %%file name /nome file
x=k(:,1);%time/tempo
y=k(:,2);%prey/prede
z=k(:,3);%predator/predatori
%Creazione della figura principale con il nome 'Lotka-Volterra'
% Create a figure with the name 'Lotka-Volterra' 
figure("Name",'Lotka-Volterra');
t1=tiledlayout(3,1) 
nexttile
%primo grafico/ first plot: time(x) and prey(y)
plot(x,y,'r') 

grid on;
xlabel('tempo');
ylabel('prede');
nexttile
%secondo grafico/ second plot: time(x) and predator(z)
plot(x,z,'g')
grid on;
xlabel('tempo');
ylabel('predatori');
nexttile
%terzo grafico/ third plot: prey(y) and predator(z)->limit cycle
plot(y,z,'c')
xlabel('prede');
ylabel('predatori');
grid on;

%Discretization/ discretizzo: second figure
figure("Name",'Lotka-Volterra quantità "naturali" ');%%variabili discretizzate
for(i=1:length(x))
    x1(i)=round(x(i)); 
    y1(i)=round(y(i));
    z1(i)=round(z(i));
end
t2=tiledlayout(3,1)
nexttile 
%primo grafico/ first plot: time(x) and prey(y)
plot(x1,y1,'r')

grid on;
xlabel('tempo');
ylabel('prede');
nexttile
%secondo grafico/ second plot: time(x) and predator(z)
plot(x1,z1,'g')
grid on;
xlabel('tempo');
ylabel('predatori');
nexttile
%terzo grafico/ third plot: prey(y) and predator(z)->limit cycle discrete
plot(y1,z1,'c')
xlabel('prede');
ylabel('predatori');
grid on;
%Log-log plot

figure("Name",'Lotka-Volterra log-log')
t3=tiledlayout(3,1)
nexttile
%primo grafico/ first plot: time(x) and prey(y)
loglog(x,y,'r');
xlabel('tempo');
ylabel('prede');
nexttile
%secondo grafico/ second plot: time(x) and predator(z)
loglog(x,z,'g')
xlabel('tempo');
ylabel('predatori');
nexttile
loglog(y,z,'c')
%terzo grafico/ third plot: prey(y) and predator(z)->limit cycle discrete
xlabel('prede');
ylabel('predatori');

%% Animation : uncomment this part
%%
    % 
    % x3=zeros(length(x),1);
    % y3=zeros(length(x),1);
    % z3=zeros(length(x),1);


% figure("Name",'Lotka-Volterra animation')
% t4=tiledlayout(3,1)
% nexttile 
% plot(x3,y3,'r')
% grid on;
% xlabel('tempo');
% ylabel('prede');
% h=animatedline;
% for(j=1:length(x))
%         x31(j)=x(j);
%         y31(j)=y(j);
%         z31(j)=z(j);
%     end
% for(j=1:length(x))
%     addpoints(h,x31(j),y31(j));
%     drawnow;
% end
% 
% 
% nexttile
% plot(x3,z3,'g');
% grid on;
% xlabel('tempo');
% ylabel('predatori');
% h2=animatedline;
% for(j=1:length(x))
%     addpoints(h2,x31(j),z31(j));
%     drawnow;
% end
% 
% nexttile
% plot(y3,z3,'c')
% xlabel('prede');
% ylabel('predatori');
% grid on;
% h2=animatedline;
% for(j=1:length(x))
%     addpoints(h2,y31(j),z31(j));
%     drawnow;
% end
