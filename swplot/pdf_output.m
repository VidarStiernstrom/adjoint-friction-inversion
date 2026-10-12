function pdf_output(f, name)
    % PDF output with embedded fonts.
    %
    %   pdf_output(f, name)          exportgraphics (vector). Fonts are embedded as
    %                                TrueType subsets. Surrounding whitespace is
    %                                trimmed, so the page is slightly smaller than
    %                                the figure (the axes keep their physical size).
    if ~endsWith(name, '.pdf')
        name = [name, '.pdf'];
    end

    use_embeddable_fonts(f)
    exportgraphics(f, name, 'ContentType', 'vector')

end

function use_embeddable_fonts(f)
    % exportgraphics leaves the 14 standard PostScript fonts unembedded (a viewer
    % is assumed to have them). Swap them for the metric-compatible system fonts,
    % which do get embedded.
    subs = { 'Helvetica',      'Arial'
             'Helvetica Neue', 'Arial'
             'Times',          'Times New Roman'
             'Courier',        'Courier New' };

    h = findall(f, '-property', 'FontName');
    for k = 1:size(subs, 1)
        hit = h(strcmpi(get(h, 'FontName'), subs{k, 1}));
        set(hit, 'FontName', subs{k, 2})
    end
end


% function pdf_output(f, name)
%     % pdf output
%     f.Units = 'centimeters';
%     f.PaperUnits = f.Units;
%     f.PaperPosition = [0, 0, f.Position(3:4)];
%     f.PaperSize = f.Position(3:4);
    
%     print(name, '-dpdf')
% end
