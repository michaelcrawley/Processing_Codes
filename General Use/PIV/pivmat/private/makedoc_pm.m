
function makedoc_pm
%MAKEDOC_PM  Generate the HTML help files for the PIVMat toolbox
%   (internal use only - 29 oct 2006)

curdir = pwd;

cd ..

makehtmldoc('*.m', 'upper', 'noad', ...
    'title', '\f (PIVMat Toolbox)',...
    'firstline', 'PIVMat Function Reference',...
    'lastline', '2005-2013 <a href="pivmat.html">PIVMat Toolbox 3.02</a>');

delete Contents.html;

movefile *.html html;

cd(curdir);

disp('HTML documentation created.');

cd(curdir);
