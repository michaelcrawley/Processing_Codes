function FixFigureFiles(flist,fig_dir)

    for n = 1:length(flist)
        filename = flist{n}(1:end-4);
        delete([fig_dir '\' filename '.png']);
        h = open([fig_dir '\' filename '.fig']);
        print(gcf,'-dpng','-r300',[fig_dir '\' filename,'.png']);
        close(h);
    end
end