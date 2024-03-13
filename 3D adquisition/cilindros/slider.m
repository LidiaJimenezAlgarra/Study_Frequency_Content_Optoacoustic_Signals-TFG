function S_f=slider(S,st,S_f,Nx,Ny,l,t,dx,dy,plano)
    if plano=='YZ'
    slider = uicontrol('Style', 'slider', 'Min', 0, 'Max', Nx, ...
        'Value', 1, 'SliderStep', [1/Ny, 1/Ny], ...
        'Position', [530 10 300 20], 'Callback', @updateImageXZ);
    else 
         slider = uicontrol('Style', 'slider', 'Min', 0, 'Max', Ny, ...
        'Value', 1, 'SliderStep', [1/Nx, 1/Nx], ...
        'Position', [530 10 300 20], 'Callback', @updateImageYZ);
    end

    function updateImageXZ(source, ~)
        value = round(get(source, 'Value')); 
        if st~='doble'
        imagesc(imrotate(squeeze(S_f(:,value,:)),-90)); colorbar; colormap('gray');
        xlabel(['Nx: ' num2str((value*dx-1e-3)*1e3) '0 mm']); title('Plano YZ');
        else
            subplot(1,2,1);imagesc(imrotate(squeeze(S(:,value,:)),-90)); colorbar; colormap('gray');
             xlabel(['Nx: ' num2str((value*dy-1e-3)*1e3) ' mm']);title('Plano YZ sin filtrar')

            subplot(1,2,2);imagesc(imrotate(squeeze(S_f(:,value,:)),-90)); colorbar; colormap('gray');
            xlabel(['Nx: ' num2str((value*dx-1e-3)*1e3) '0 mm']); title('Plano YZ filtrado');
        end
    end
    function updateImageYZ(source, ~)
        value = round(get(source, 'Value')); % Obtiene el valor del slider
        if st~='doble'
          imagesc(imrotate(squeeze(S_f(value,:,:)),-90)); colorbar; colormap('gray');
          xlabel(['Ny: ' num2str((value*dx-2e-3)*1e3) '0 mm']); title('Plano XZ');
        else
        subplot(1,2,1);imagesc(imrotate(squeeze(S(value,:,:)),-90)); colorbar; colormap('gray');
        xlabel(['Ny: ' num2str((value*dx-2e-3)*1e3) '0 mm']); title('Plano XZ sin filtrar');

        subplot(1,2,2);imagesc(imrotate(squeeze(S_f(value,:,:)),-90)); colorbar; colormap('gray');
        xlabel(['Ny: ' num2str((value*dx-2e-3)*1e3) '0 mm']); title('Plano XZ filtrado');
        end
    end
end
     