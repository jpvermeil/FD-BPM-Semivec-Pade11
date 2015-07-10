function [Px,glob_adr_slgs,dim_xl,dim_yl] = FD_BPM_Pade11_Semivec_vec_(n,lambda,n_eff_Px,n_eff_Py,alpha,solver_tol,xg,yg,dz,Nb_interpolation,POLARIZATION,FIELDCOMPONENTS,BC)

% BPM zur Bestimmung von E- oder H-Feldern. Bestimmung erfolgt 
% komponentenweise. Die Polarisation kann mit POLARIZATION (String) 
% ausgewählt werden und kann zu 'TE' oder 'TM' gesetzt werden.
% Nb_interpolation setzt die ver-2^Nb_interpolationsfacht das
% Brechzahlprofil in beide Richtungen vor der Bestimmung der
% Ausbreitungskonstanten. 

%% Wide angle Pade(3,3) FD BPM for E- and H fields

format long

%% Variablen und Konstanten definieren

beta_0 = 2*pi/lambda;
beta_z = beta_0 * n_eff_Px;

% Feldkomponenten werden überprüft und festgelegt

[FX,FY] = check_fieldcomponents_(FIELDCOMPONENTS);
        
%% Hilfsgrößen

dim_y   = size(n,1);
dim_x   = size(n,2);

dim_xl  = size(n,2) - 2;
dim_yl  = size(n,1) - 2;

Px = zeros(size(n,1),size(n,2),size(n,3));
Py = zeros(size(n,1),size(n,2),size(n,3));

dG.cl = zeros(size(n,1),size(n,2));
dG.cg = zeros(size(n,1),size(n,2));

dG.f = zeros(size(n,1),size(n,2));
dG.f([1 end],[1:end 1:end]) = 1;
dG.f([1:end 1:end],[1 end]) = 1;

dG.cl(2:end-1,2:end-1) = reshape(1:1:dim_xl*dim_yl',dim_yl,dim_xl);
dG.cg(1:end) = 1:1:length(dG.cg(1:end));

glob_adr_slgs   = dG.cg(2:end-1,2:end-1);
glob_adr_slgs   = reshape(glob_adr_slgs,size(glob_adr_slgs,1)*size(glob_adr_slgs,2),1);
     
%% Interpolation des Brechzahlprofils zur Bestimmung der Ausbreitungskonstanten und des Grundmodenfeldes

ni = interpn(squeeze(n(:,:,1)),Nb_interpolation);    % ver-Nb_interpolation-fachung des Brechzahlprofils

%% Bestimmung der Ausbreitungskonstanten und des Grundmodenfeldes

[~,Modenfeld_Px] = FD_propagationconstants_Semivec_(ni,beta_0,xg(:,:,1),yg(:,:,1),dim_y,dim_xl,dim_yl,dG.cg,dG.cl,POLARIZATION,FX,1); 
%[~,Modenfeld_Py] = FD_propagationconstants_Semivec_(ni,beta_0,xg(:,:,1),yg(:,:,1),dim_y,dim_xl,dim_yl,dG.cg,dG.cl,POLARIZATION,FY,1);

Px_a = zeros(size(ni,1),size(ni,2));
%Py_a = zeros(size(ni,1),size(ni,2));

Px_a(2:end-1,2:end-1) = reshape(Modenfeld_Px(:,1),dim_yl,dim_xl);
%Py_a(2:end-1,2:end-1) = reshape(Modenfeld_Py(:,1),dim_yl,dim_xl);

Px_a = Px_a/max(max(abs((Px_a))));
%Py_a = Py_a/max(max(abs((Py_a))));

%% Anfangswerte für die Felder

Px(:,:,1) = Px_a;

xg1 = squeeze(xg(:,:,1));
yg1 = squeeze(yg(:,:,1));

[r_max,c_max] = find(squeeze(n(:,:,1)) == max(max(squeeze(n(:,:,1)))));

xg1 = xg1 - xg1(r_max(1),c_max(1));
yg1 = yg1 - yg1(r_max(1),c_max(1));

test = 1*exp(-xg1.^2/(3e-6)^2 -yg1.^2/(3e-6).^2);
[r_max,c_max] = find(test == max(max(test)));

% figure
% surf(xg1,yg1,test)
% shading flat

Px_a(:,:) = test;
Px(:,:,1) = Px_a;

% n_threshold = min(min(squeeze(n(:,:,1)))) + (.5)*(max(max(squeeze(n(:,:,1)))) - min(min(squeeze(n(:,:,1)))));
% Px_a = zeros(size(ni,1),size(ni,2));
% Px_a(find(squeeze(n(:,:,1)) >= n_threshold)) = 1;
% %Px_a(:,:) = test;
% Px(:,:,1) = Px_a;
% surf(xg(:,:,1),yg(:,:,1),Px_a)

%% Multistep Variablen generieren

% Die Multistep Variablen sind für gleich bleibende Propagationsschritte
% und effektive Brechzahlen über die gesamte Propagationsstrecke gleich und
% müssen nicht erneut berechnet werden. Ändert sich einer dieser beiden
% Parameter, bei zum Beispiel einem nicht äquidistanten Propagationsgrid,
% müssen diese in der folgende Schleife mitberechnet werden.

[ux,vx] = gen_multistep_vars_11_(dz,alpha,beta_z);   
%[uy,vy] = gen_multistep_vars_(dz,alpha,n_eff_Py);

%% Propagation in z-Richtung

tic
c_alt = 1;
c = 1; % Globaler Forschrittszähler

h = waitbar(0,'','Name',['Berechne Padé ' POLARIZATION '-BPM mit ' BC '-Randbedingung...']);

for kz = 1:1:size(n,3)-1
    
    %% Rechte Seite erstellen
    
    % Extrahiere bekanntes Feld bei kz in zweidimensionale Matrix
    
    Pbx = Px(:,:,kz);
    
    %% Randbegindung
    
    if strcmp(BC,'TBC')
       
        % In diesem Schritt werden die bekannten Feldwerte der Berandung so
        % angepasst, dass Leistung reflexionsfrei aus dem Gebiet austreten
        % können ohne reflektiert zu werden. Dies ist insbesondere deshalb
        % wichtig, weil der Großteil der Leistung bei diffundierten
        % Wellenleiterstrukturen nah an der oberen Berandung geführt wird. Je
        % tiefer der Wellenleiter vergragen ist, desto weniger kritisch ist die
        % Betrachtung der insgesamt geführten Leistung. Da Dämpfungsverhalten
        % untersucht wird würden also reflektierte Anteile das Ergebnis
        % verfälschen.

        % Upper boundary

        kx = 2:1:size(n,2)-1;
        dy_u = abs(yg(end-1,kx,kz) - yg(end-2,kx,kz));
        eta_u = Pbx(end-1,kx)./Pbx(end-2,kx);
        gamma_tbc{1} = eta_u;
        beta_x_u = log(eta_u)./(-1j*dy_u);    
        beta_x_u(find(real(beta_x_u) > 0)) = beta_x_u(find(real(beta_x_u) > 0)) - 2*real(beta_x_u(find(real(beta_x_u) > 0)));
        %Pbx(end-1,kx) = Pbx(end-2,kx) .* exp(-1j*beta_x_u.*dy_u);
        %Pbx(end,kx) = Pbx(end-1,kx) .* exp(-1j*beta_x_u.*dy_u);
        
        eta_u = Pbx(end,kx)./Pbx(end-1,kx);
        gamma_tbc{1} = eta_u;
        beta_x_u_c = log(eta_u)./(-1j*dy_u);    
        
        if isempty(real(beta_x_u_c) > 0)
            
            out = 'WARNING: Positive elements in transversal wavenumber beta_x_u after correcting hypothetical boundary nodes!';
            disp(out);
            return
        
        end
        
        % Right hand boundary

        ky = 2:1:size(n,1)-1;
        dx_r = abs(xg(ky,end-1,kz) - xg(ky,end-2,kz));
        eta_r = Pbx(ky,end-1)./Pbx(ky,end-2);
        gamma_tbc{3} = eta_r;
        beta_x_r = log(eta_r)./(-1j*dx_r);    
        beta_x_r(find(real(beta_x_r) > 0)) = beta_x_r(find(real(beta_x_r) > 0)) - 2*real(beta_x_r(find(real(beta_x_r) > 0)));
        %Pbx(ky,end-1) = Pbx(ky,end-2) .* exp(-1j*beta_x_r.*dx_r);
        %Pbx(ky,end) = Pbx(ky,end-1) .* exp(-1j*beta_x_r.*dx_r);

        eta_r = Pbx(ky,end)./Pbx(ky,end-1);
        gamma_tbc{3} = eta_r;
        beta_x_r_c = log(eta_r)./(-1j*dx_r); 
        
        if isempty(real(beta_x_r_c) > 0)
            
            out = 'WARNING: Positive elements in transversal wavenumber beta_x_r after correcting hypothetical boundary nodes!';
            disp(out);
            return
        
        end
        
        % Lower boundary

        kx = 2:1:size(n,2)-1;
        dy_l = abs(yg(3,kx,kz) - yg(2,kx,kz));
        eta_d = Pbx(3,kx)./Pbx(2,kx);
        gamma_tbc{2} = 1./eta_d;
        beta_x_d = log(eta_d)./(1j*dy_l);    
        beta_x_d(find(real(beta_x_d) > 0)) = beta_x_d(find(real(beta_x_d) > 0)) - 2*real(beta_x_d(find(real(beta_x_d) > 0)));
        %Pbx(2,kx) = Pbx(3,kx) .* exp(-1j*beta_x_d.*dy_l);
        %Pbx(1,kx) = Pbx(2,kx) .* exp(-1j*beta_x_d.*dy_l);
        
        eta_d = Pbx(2,kx)./Pbx(1,kx);
        gamma_tbc{2} = 1./eta_d;
        beta_x_d_c = log(eta_d)./(1j*dy_l);   
        
        if isempty(real(beta_x_d_c) > 0)
            
            out = 'WARNING: Positive elements in transversal wavenumber beta_x_d after correcting hypothetical boundary nodes!';
            disp(out);
            return
        
        end
        
        % Left hand boundary

        ky = 2:1:size(n,1)-1;                                                                                                       % Vektor der zu untersuchenden Berandungselemente
        dx_l = abs(xg(ky,3,kz) - xg(ky,2,kz));
        eta_l = Pbx(ky,3)./Pbx(ky,2);                                                                                         % Verhältnis der bekannten Felder
        gamma_tbc{4} = 1./eta_l;                                                                                                         % Verhältnis der bekannten Felder (Koeffizientendarstellung)
        beta_x_l = log(eta_l)./(1j*dx_l);                                                                                            % Transversale Wellenvektorkomponente
        beta_x_l(find(real(beta_x_l) > 0)) = beta_x_l(find(real(beta_x_l) > 0)) - 2*real(beta_x_l(find(real(beta_x_l) > 0)));       % Anpassung des Realteils (muss negativ sein, damit Leistung reflexionsfrei austreten kann!)
        %Pbx(ky,2) = Pbx(ky,3) .* exp(-1j*beta_x_l.*dx_l);
        %Pbx(ky,1) = Pbx(ky,2) .* exp(-1j*beta_x_l.*dx_l); % Gl. 5.166   % Neuberechnung des bekannten Feldwertes auf der Berandung
        
        eta_l = Pbx(ky,2)./Pbx(ky,1);                                                                                         % Verhältnis der bekannten Felder
        gamma_tbc{4} = 1./eta_l;  
        beta_x_l_c = log(eta_l)./(1j*dx_l); 
        
        if isempty(real(beta_x_l_c) > 0)
            
            out = 'WARNING: Positive elements in transversal wavenumber beta_x_l after correcting hypothetical boundary nodes!';
            disp(out);
            return
        
        end
        
    elseif strcmp(BC,'ABC')
        
        gamma_tbc = 0;
    
    end
    
    %% Assemblieren des LGS mit Vektoransatz (etwa 5 Mal schneller als mit Schleifen!)
     
    % Es ist an dieser Stelle anzumerken, dass das Erstellen der
    % Sparse-Matrix mit for Schleifen und auch mit dem sparse Befehl zu
    % lange dauert um Berechnungen dieser Größe in sinnvollen Zeiten
    % durchführen zu können. spdiags eignet sich am besten für die
    % Zuweisung der per Funktionsaufruf generierten Diagonalen. Der Gewinn
    % von spdiags gegenüber sparse liegt bei Faktor 50 bzw. 250 im
    % Vergleich zur doppelten for-schleife. spdiags ist die einzige
    % Variante die ernsthaft zur Berechnung benutzt werden sollte.
    % Das erstellen der Systemmatrix als Sparse Matrix ist die
    % zeitaufwändigste Operation jedes Propagationsschrittes (für dim_slgs < 10k x 10k). Vorherige
    % Allokation bringt keinen zeitlichen Gewinn gegenüber der direkten
    % Zuweisung der Diagonalen. Für große Gleichungssysteme (dim_slgs > 10k
    % x 10k) dauert die eigentliche Lösung wesentlich länger als die
    % Erstellung der Systemmatrix.
 
%% Brechzahlprofil für ersten ADI-Schritt:

    %n_adi_1 = (n(:,:,kz+1) + n(:,:,kz))/2;
    
%% Funktionen zur Diagonalenerstellung ausführen (nur einmal pro kz nötig)
    
    [ Cxx,Nxx,Sxx,Exx,Wxx ]      = diagonals_pade_(beta_0,n_eff_Px,n(:,:,kz+1),xg(:,:,kz+1),yg(:,:,kz+1),dim_y,dim_xl,dim_yl,dG,gamma_tbc,POLARIZATION,FX,BC);
    [ Cbxx,Nbxx,Sbxx,Ebxx,Wbxx ] = diagonals_pade_(beta_0,n_eff_Px,n(:,:,kz),xg(:,:,kz),yg(:,:,kz),dim_y,dim_xl,dim_yl,dG,gamma_tbc,POLARIZATION,FX,BC);
    
    %% Multistep Methode anwenden. 3 Schritte für Padé(3,3) nötig
    
    for ii = 1:1:1
        
        % Linke Seite bleibt bis auf den Multistepparameter u(ii) gleich für kz = const. 
        
        % x-Komponente des Feldes
        
        Axx = sparse((size(n,1)-2)*(size(n,2)-2),(size(n,1)-2)*(size(n,2)-2));
        Axx = spdiags(1 + vx(ii)*Cxx,0,Axx);
        Axx = spdiags([vx(ii)*Nxx(2:end); 0],-1,Axx);
        Axx = spdiags([0; vx(ii)*Sxx(1:end-1)],1,Axx);
        Axx = spdiags([zeros(dim_yl,1); vx(ii)*Exx(1:end-dim_yl)],dim_yl,Axx);
        Axx = spdiags([vx(ii)*Wxx(dim_yl+1:end); zeros(dim_yl,1)],-dim_yl,Axx);
        
        Cbxx = 1 + ux(ii)*Cbxx;
        Nbxx = [ux(ii)*Nbxx(2:end); 0];
        Sbxx = [0; ux(ii)*Sbxx(1:end-1)];
        Ebxx = [zeros(dim_yl,1); ux(ii)*Ebxx(1:end-dim_yl)];
        Wbxx = [ux(ii)*Wbxx(dim_yl+1:end); zeros(dim_yl,1)];
        
        Cbxx = Pbx(glob_adr_slgs)          .* Cbxx;
        Nbxx = Pbx(glob_adr_slgs - 1)      .* Nbxx;
        Sbxx = Pbx(glob_adr_slgs + 1)      .* Sbxx;
        Ebxx = Pbx(glob_adr_slgs + dim_y)  .* Ebxx;
        Wbxx = Pbx(glob_adr_slgs - dim_y)  .* Wbxx;
        
        bxx = sparse(Cbxx + Nbxx + Sbxx + Ebxx + Wbxx);

        

%         Bxx = sparse((size(n,1)-2)*(size(n,2)-2),(size(n,1)-2)*(size(n,2)-2));
%         Bxx = spdiags(1 + vx(ii)*Cbxx,0,Bxx);
%         Bxx = spdiags([vx(ii)*Nbxx(2:end); 0],-1,Bxx);
%         Bxx = spdiags([0; vx(ii)*Sbxx(1:end-1)],1,Bxx);
%         Bxx = spdiags([zeros(dim_yl,1); vx(ii)*Ebxx(1:end-dim_yl)],dim_yl,Bxx);
%         Bxx = spdiags([vx(ii)*Wbxx(dim_yl+1:end); zeros(dim_yl,1)],-dim_yl,Bxx);
%         
%         Pbxx = reshape(Pbx(2:end-1,2:end-1),1,dim_xl*dim_yl)';
% 
%         bxx = Bxx * Pbxx;

%         Bxx = sparse((size(n,1)-2)*(size(n,2)-2),(size(n,1)-2)*(size(n,2)-2));
%         Bxx = spdiags(1 + ux(ii)*Cxx,0,Bxx);
%         Bxx = spdiags([ux(ii)*Nxx(2:end); 0],-1,Bxx);
%         Bxx = spdiags([0; ux(ii)*Sxx(1:end-1)],1,Bxx);
%         Bxx = spdiags([zeros(dim_yl,1); ux(ii)*Exx(1:end-dim_yl)],dim_yl,Bxx);
%         Bxx = spdiags([ux(ii)*Wxx(dim_yl+1:end); zeros(dim_yl,1)],-dim_yl,Bxx);
%         
%         bxx = Bxx * sparse(reshape(Pbx(2:end-1,2:end-1),1,dim_xl*dim_yl)');
        
%         Cbxx = 1 + ux(ii)*Cxx;
%         Nbxx = [ux(ii)*Nxx(2:end); 0];
%         Sbxx = [0; ux(ii)*Sxx(1:end-1)];
%         Ebxx = [zeros(dim_yl,1); ux(ii)*Exx(1:end-dim_yl)];
%         Wbxx = [ux(ii)*Wxx(dim_yl+1:end); zeros(dim_yl,1)];
%         
%         bx = sparse(Cbxx + Nbxx + Sbxx + Ebxx + Wbxx);
%         bx = sparse(reshape(Pbx(2:end-1,2:end-1),1,dim_xl*dim_yl)' .* bx);
        
        % y-Komponente des Feldes
        
%         Ayy = sparse((size(n,1)-2)*(size(n,2)-2),(size(n,1)-2)*(size(n,2)-2));
%         Ayy = spdiags(1 + vy(ii)*Cyy,0,Ayy);
%         Ayy = spdiags([vy(ii)*Nyy(2:end); 0],-1,Ayy);
%         Ayy = spdiags([0; vy(ii)*Syy(1:end-1)],1,Ayy);
%         Ayy = spdiags([zeros(dim_yl,1); vy(ii)*Eyy(1:end-dim_yl)],dim_yl,Ayy);
%         Ayy = spdiags([vy(ii)*Wyy(dim_yl+1:end); zeros(dim_yl,1)],-dim_yl,Ayy);
%         
%         Cbyy = 1 + uy(ii)*Cyy;
%         Nbyy = [uy(ii)*Nyy(2:end); 0];
%         Sbyy = [0; uy(ii)*Syy(1:end-1)];
%         Ebyy = [zeros(dim_yl,1); uy(ii)*Eyy(1:end-dim_yl)];
%         Wbyy = [uy(ii)*Wyy(dim_yl+1:end); zeros(dim_yl,1)];
%         
%         by = sparse(Cbyy + Nbyy + Sbyy + Ebyy + Wbyy);
%         by = sparse(reshape(Pby(2:end-1,2:end-1),1,dim_xl*dim_yl)' .* by);
        

        %% Lösung des Gleichungssystems
        
        [Px_l,~] = bicgstab(Axx,bxx,solver_tol);
        %Px_l = Axx\bxx;
        %[Py_l,~] = bicgstab(Ayy,by,solver_tol);
        
        Pbx = zeros(dim_y,dim_x); 
        Pbx(2:end-1,2:end-1) = reshape(Px_l,dim_yl,dim_xl);
        
        %% Absorber hinzufügen
        
         adr_n_threshold         = find(squeeze(n(:,:,kz)) <= 1.6);
         Pbx(adr_n_threshold)    = 0;

    end
    
    %% Rekonstruktion der Lösung
    
    Px(:,:,kz+1) = Pbx;
    
%     s = 5;                                                      % Prozentinkrement zur Fortschrittsmeldung. s = 5 -> Fortschrittsmeldung alle 5%.
%     if mod(kz,round(size(n,3)/(100/s))) == 0 
%         
%         Dauer_Minuten     = floor(toc/60);
%         Dauer_Sekunden    = toc - Dauer_Minuten*60;
%         Verb_Minuten      = floor((toc/c_alt)*((100/s)-c_alt)/60);
%         Verb_Sekunden     = (toc/c_alt)*((100/s)-c_alt) - Verb_Minuten*60;
%         out = ['Fortschritt: ' num2str(c_alt*s) '% nach bisher ' num2str(Dauer_Minuten) 'm ' num2str(Dauer_Sekunden) 's. Verbleibende Dauer ca. ' num2str(Verb_Minuten) 'm ' num2str(Verb_Sekunden) 's.'];
%         disp(out);
%         c_alt = c_alt + 1;
%         
%     end
    
    
%     if ~mod(kz,round(size(n,3)/100)) 
%       % 
%       % Es kann passieren, das der Benutzer die WAITBAR mit der Mouse-Click 
%       % auf den X-Button abschiesst. Um sicher zu sein, dass die WAITBAR 
%       % noch da ist, testen wir, on der die noch aktiv ist 
%       if ishandle(h) 
%          % 
%          % die eigentliche WAITBAR Aktualisierung. 
%          % Nebenbei wird auch der erreichter Stand in "%" dargestellt. 
%          
%         Dauer_Minuten     = floor(toc/60);
%         Dauer_Sekunden    = toc - Dauer_Minuten*60;
%         Verb_Minuten      = floor((toc/c)*((100/s)-c)/60);
%         Verb_Sekunden     = (toc/c)*((100/s)-c) - Verb_Minuten*60;
%         %out = ['Fortschritt: ' num2str(c*s) '% nach bisher ' num2str(Dauer_Minuten) 'm ' num2str(Dauer_Sekunden) 's. Verbleibende Dauer ca. ' num2str(Verb_Minuten) 'm ' num2str(Verb_Sekunden) 's.'];
%         %disp(out);
%         c = c + 1;
%         
%         waitbar(kz/size(n,3),h,['Fortschritt: ' num2str(c*s) '%. Etwa verbleibend: ' num2str(Verb_Minuten) 'm ' num2str(Verb_Sekunden) 's.']) 
%       
%       else 
%          % 
%          disp('Prozess wurde unterbrochen') 
%          break 
%       end 
%    end 
    s = 1;
    if floor(100*kz/(size(n,3)-1)) > c*s 
      % 
      % Es kann passieren, das der Benutzer die WAITBAR mit der Mouse-Click 
      % auf den X-Button abschiesst. Um sicher zu sein, dass die WAITBAR 
      % noch da ist, testen wir, on der die noch aktiv ist 
      if ishandle(h) 
         % 
         % die eigentliche WAITBAR Aktualisierung. 
         % Nebenbei wird auch der erreichter Stand in "%" dargestellt. 
         
        Dauer_Minuten     = floor(toc/60);
        Dauer_Sekunden    = toc - Dauer_Minuten*60;
        Verb_Minuten      = floor((toc/c)*((100/s)-c)/60);
        Verb_Sekunden     = (toc/c)*((100/s)-c) - Verb_Minuten*60;
        %out = ['Fortschritt: ' num2str(c*s) '% nach bisher ' num2str(Dauer_Minuten) 'm ' num2str(Dauer_Sekunden) 's. Verbleibende Dauer ca. ' num2str(Verb_Minuten) 'm ' num2str(Verb_Sekunden) 's.'];
        %disp(out);
        c = c + 1;
        
        waitbar(kz/(size(n,3)-1),h,['Fortschritt: ' num2str(c*s) '%. Etwa verbleibend: ' num2str(Verb_Minuten) 'm ' num2str(Verb_Sekunden) 's.']) 
      
      else 
         % 
         disp('Prozess wurde unterbrochen') 
         break 
      end 
   end 
    
    
    
end

try
    close(h)
end

Dauer_Minuten     = floor(toc/60);
Dauer_Sekunden    = toc - Dauer_Minuten*60;

out = ['BPM-Berechnung beendet. Gesamte Dauer: ' num2str(Dauer_Minuten) 'm ' num2str(Dauer_Sekunden) 's.'];
disp(out);

end
