function [rgb hex] = wes_palettes(palettename)
%colors = wes_palettes(palettename)
%get color palettes from Wes Anderson movies. Available:
% 'GrandBudapest', 'Moonrise', 'Royal', 'Moonrise2', 'Cavalcanti',
% 'Royal2', 'GrandBudapest2', 'Moonrise3', 'Chevalier', 'Zissou',
% 'FantasticFox', 'Darjeeling', 'Rushmore', 'BottleRocket', 'Darjeeling2',
% 'BottleRocket2', 'Royal3', 'Royal4', 'GrandBudapest3', 'GrandBudapest4',
% 'GrandBudapest5', 'Zissou2', 'Zissou3', 'Moonrise4', 'Moonrise5',
% 'FantasticFox2', 'Darjeeling3'

if nargin < 1
    palettename = 'GrandBudapest';
    disp('Usage: colors = wes_palettes(palettename);');
    disp('get color palettes from Wes Anderson movies. Available:');
    disp(' GrandBudapest, Moonrise, Royal, Moonrise2, Cavalcanti,');
    disp(' Royal2, GrandBudapest2, Moonrise3, Chevalier, Zissou,');
    disp(' FantasticFox, Darjeeling, Rushmore, BottleRocket, Darjeeling2,');
    disp(' BottleRocket2, Royal3, Royal4, GrandBudapest3, GrandBudapest4,');
    disp(' GrandBudapest5, Zissou2, Zissou3, Moonrise4, Moonrise5,');
    disp(' FantasticFox2, Darjeeling3');
end

%assign color palettes (hex values)
% these are from an R script online
wes = struct;
wes.GrandBudapest  = {'#F1BB7B', '#FD6467', '#5B1A18', '#D67236'};
wes.GrandBudapest2 = {'#E6A0C4', '#D8A499', '#C6CDF7', '#7294D4'};

wes.Moonrise       = {'#F3DF6C', '#CEAB07', '#D5D5D3', '#24281A'};
wes.Moonrise2      = {'#798E87', '#C27D38', '#CCC591', '#29211F'};
wes.Moonrise3      = {'#85D4E3', '#F4B5BD', '#9C964A', '#CDC08C', '#FAD77B'};

wes.Royal          = {'#899DA4', '#C93312', '#FAEFD1', '#DC863B'};
wes.Royal2         = {'#9A8822', '#F5CDB4', '#F8AFA8', '#FDDDA0', '#74A089'};

wes.Darjeeling     = {'#FF0000', '#00A08A', '#F2AD00', '#F98400', '#5BBCD6'};
wes.Darjeeling2    = {'#ECCBAE', '#046C9A', '#D69C4E', '#ABDDDE', '#000000'};

wes.BottleRocket   = {'#A42820', '#5F5647', '#9B110E', '#3F5151', '#4E2A1E', '#550307', '#0C1707'};
wes.BottleRocket2  = {'#FAD510', '#CB2314', '#273046', '#354823', '#1E1E1E'};

wes.Cavalcanti     = {'#D8B70A', '#02401B', '#A2A475', '#81A88D', '#972D15'};
wes.Chevalier      = {'#446455', '#FDD262', '#D3DDDC', '#C7B19C'};
wes.Zissou         = {'#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00'};
wes.FantasticFox   = {'#DD8D29', '#E2D200', '#46ACC8', '#E58601', '#B40F20'};
wes.Rushmore       = {'#E1BD6D', '#EABE94', '#0B775E', '#35274A' ,'#F2300F'};

% these are from wesandersonpalettes.tumblr.com
wes.Royal3         = {'#FE6D7E', '#F19167', '#E8AA9D', '#3D1F1F', '#FCDED4'};
wes.Royal4         = {'#99B1B3', '#D3AF7B', '#DED4CB', '#DEB0B0'};

wes.GrandBudapest3 = {'#E7A396', '#9F295D', '#D9A34B', '#ECD1CA', '#430E00'};
wes.GrandBudapest4 = {'#C68006', '#642E68', '#FFF7C9', '#A9902B', '#2A1824'};
wes.GrandBudapest5 = {'#FFE2C4', '#314856', '#ECBBAD', '#798881', '#9D5141', '#A67F48'};

wes.Zissou2        = {'#9DADC4', '#EBE85D', '#6CAFCC', '#AF9E73', '#C9C6D7'};
wes.Zissou3        = {'#433022', '#235135', '#F3D28F', '#7B6DA8', '#C5A495'};

wes.Moonrise4      = {'#CB654F', '#D3B1A7', '#CFCB9C', '#8CBEA3', '#DFBA47'};
wes.Moonrise5      = {'#E69EA1', '#DED9A2', '#F6C83E', '#4D5B28', '#DB4372', '#B77E60'};

wes.FantasticFox2  = {'#E8C95D', '#D16B54', '#A9D8C8', '#433447', '#B9B09F'};
wes.Darjeeling3    = {'#DDE8E0', '#749CA8', '#F9E0A8', '#BC5E21', '#8B8378'};

%check if movie is available
fnames = fieldnames(wes);
if ~any(strcmp(fnames,palettename))
    palettes = fieldnames(wes);
    fprintf('\nAvailable palettes:\n');
    for ip = 1:length(palettes)
        fprintf('\n%s',palettes{ip});
    end
    fprintf('\n\n');
    error('wes_palettes: Unknown Palette Name! Choose one of the above.')
end

%get hex values for selected movie
hex = wes.(palettename);
if isrow(hex), hex = hex'; end
hex = cell2mat(hex);
if strcmpi(hex(1,1),'#'), hex(:,1) = []; end

%convert hex to rgb
rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;