function make_movie(paths, plot_variables, uvmaxs, Nres, sim_types, axlims, fault_xlims)
    default_arg('plot_variables', 'all');
    default_arg('uvmaxs', struct());
    default_field(uvmaxs, 'true', [1, 1.5]);
    default_field(uvmaxs, 'fwd', [1, 1.5]);
    default_field(uvmaxs, 'adj', [0.2, 0.2]);
    default_arg('Nres', 69); % number of results to plot
    default_arg('sim_types', {'true', 'fwd', 'adj'});
    default_arg('axlims', [-10, 10, -10, 10]); % xleft, xright, ybottom, ytop
    default_arg('fault_xlims', [-10, 10]); % xleft, xright
    
    type = struct();
    type.plot_variables = plot_variables; % 'all' or 'vVS'
    % type.axlims = [-10, 10, -10, 10]; % xleft, xright, ybottom, ytop
    type.axlims = axlims;
    type.fault_xlims = fault_xlims;

    for i = 1:length(sim_types)
        sim_type = sim_types{i};
        % savedir =  sprintf('%s%s/%s', outpath, name, sim_type);
        savedir = paths.(sim_type);
        ini_path = sprintf('%s/%d.mat',savedir,1);

        type.umax = uvmaxs.(sim_type)(1);
        type.vmax = uvmaxs.(sim_type)(2);
        if sim_type == "adj"
            type.sim = "adj";
        else
            type.sim = "fwd";
        end 
        obj = PlotSW(ini_path, type);

        % vname = sprintf('%s%s/%s.mp4', outpath, name, sim_type);
        if ispc || ismac
            vname = sprintf('%s/../%s.mp4', savedir, sim_type);
            mov = VideoWriter(vname, 'MPEG-4');
        else % MPEG-4 is not supported by VideoWriter on Linux
            vname = sprintf('%s/../%s.avi', savedir, sim_type);
            mov = VideoWriter(vname, 'Motion JPEG AVI');
        end
        mov.FrameRate = 5;
        open(mov);
        for it = 1:Nres
            fpath = sprintf('%s/%d.mat',savedir,it);
            obj.plot_result(fpath);
            set(gcf,'Color','w');
            frame = getframe(gcf);
            writeVideo(mov, frame);
        end
        close(mov);
    end
end