function VisuIdB(varargin)
% VisuIdB:  visualise "l'intensité" en dB d'une image en amplitude complexe sur une échelle de gris (English help below)
% USAGE: VisuIdB(ImAC,seuil_dB)  ou   VisuIdB(ImAC,seuil_dB,tronc_dB)  ou
%        VisuIdB(x,y,ImAC,seuil_dB)  ou   VisuIdB(x,y,ImAC,seuil_dB,tronc_dB)
% x       : [OPTIONNEL] vecteur de 2 (ou n) éléments caractérisant les bornes horizontales de l'image à afficher
% y       : [OPTIONNEL] vecteur de 2 (ou m) éléments caractérisant les bornes verticales de l'image à afficher
%           Si x a plus que 2 éléments, seuls le premier et dernier élément est utilisé x=[x(1) x(end)]. Idem pour y.
%           ## x et y doivent être SIMULTANÉMENT présents ou absents. En cas d'absence, les valeurs par DÉFAUT
%           sont x=[1 n] et y=[1 m], où [m,n]=size(ImAC).
% ImAC    : matrice complexe (taille m×n) contenant l'image en amplitude complexe à afficher
% seuil_dB: scalaire réelle indiquant le seuil bas en dB de l'image affichée (typiquement -30 ou -40). Toutes valeurs
%           de l'intensité inférieures ou égales à ce seuil sera affichées en noir. (Le maximum est affiché en blanc)
% tronc_dB: [OPTIONNEL]. Si absent, le seuil précédent "seuil_dB" est RELATIF AU MAXIMUM de l'image.
%                        Le max de l'image est normalisé à 0dB pour l'affichage.
%                        Si présent, le seuil précédent "seuil_dB" est ABSOLU et "tronc_dB" indique un niveau
%                        de troncature haute au-delà duquel l'affichage est blanc (Indiquer +Inf pour supprimer cette
%                        troncature et juste afficher des dB absolus et non pas relatifs aux maximum).


switch nargin
  
  case {2,3} % cas VisuIdB(ImAC,seuil_dB) ou VisuIdB(ImAC,seuil_dB,tronc_dB)
    AC=varargin{1};
    seuil_dB=varargin{2};
    if nargin==3
      tronc_dB=varargin{3};
      if ~isscalar(tronc_dB)
        error(['La syntaxe "VisuIdB(x,y,ImAC)" est invalide. Il faut préciser "seuil_dB"' char(10) ...
               'The syntax "VisuIdB(x,y,ImAC)" is incorrect. A fourth argument "Threshold_dB" is required']) %#ok<*CHARTEN>
      end
      tronc=true;
    else
      tronc=false;
    end
    [m,n]=size(AC);
    x=[1 n];
    y=[1 m];
    YDir='reverse';
    
  case {4,5} % cas VisuIdB(x,y,ImAC,seuil_dB)  ou   VisuIdB(x,y,ImAC,seuil_dB,tronc_dB)
    x=varargin{1};
    y=varargin{2};
    YDir='normal';
    AC=varargin{3};
    seuil_dB=varargin{4};
    if nargin==5
      tronc_dB=varargin{5};
      tronc=true;
    else
      tronc=false;
    end
    
  otherwise
    error(['Syntaxe d''appel incorrect! Il faut entre 2 et 5 arguments en entrée; voir l''aide-en-ligne.' char(10) ...
           'Incorrect calling syntax! Two to five input arguments are required; see on-line help.'])
    
end

U=real(AC.*conj(AC));

seuil_eff=10^(seuil_dB/10);

if ~tronc
  % seuil_dB RELATIF:
  if seuil_dB>=0
    error(['L''argument "seuil_dB" doit être négatif (valeur relative au max à 0dB)' char(10) ...
           'The "Thres_dB" argument must be negative (value relative to the max set at 0dB)'])
  end
  figure(gcf),imagesc(x,y,10*log10(max(seuil_eff,U/max(U(:))))),colormap(gca,gray(256)),axis image,h=colorbar;
  % log pris après le seuillage pour éviter un éventuel warning log of 0.
  set(get(h,'Title'),'String','dB_{relatif}')
  set(get(h,'XLabel'),'String',sprintf('max=%1.1e',max(U(:))),'FontSize',8)
else
  figure(gcf),imagesc(x,y,min(tronc_dB,10*log10(max(seuil_eff,U)))),colormap(gca,gray(256)),axis image,h=colorbar;
  % log pris après le seuillage pour éviter un éventuel warning log of 0.
  set(get(h,'Title'),'String','dB_{absolu}')
end
set(gca,'YDir',YDir)
set(gca,'TickDir','out')
set(h,'TickDir','out')
