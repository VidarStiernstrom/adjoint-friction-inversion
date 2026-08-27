function pdf_output(f, name)
    % pdf output
    f.Units = 'centimeters';
    f.PaperUnits = f.Units;
    f.PaperPosition = [0, 0, f.Position(3:4)];
    f.PaperSize = f.Position(3:4);
    
    print(name, '-dpdf')
end
