function md = initParams(setup, varargin )
% Initializes metadata structure with expected set of parameters.


    ref = eval(setup);

    %Copy requested reference parameters to md

    if nargin>1
        if isstruct( varargin{1} )
            % If the first input is a structure, append parameters
            md = varargin{1};
            first = 2;
        else
            % Otherwise start from scratch
            md = struct();
            first = 1;
        end

        for i=first:nargin
            var = varargin{i};

            if isfield( md, var )
                warning( 'Paramter already initialized: %s', var );
            elseif isfield( ref, var )
                md.(var) = ref.(var);
            else
                md.(var) = parameter;
                warning( 'Unrecognized parameter: %s', var );
            end
        end
    else
        % If no paramters are requested, output all of them
        md = ref;
    end

end

function ref = NearFieldParams()
    %Reference parameters   
    ref.M			= parameter( 'Mach Number', 'M', '' );
    ref.m           = parameter( 'Forcing Mode','m', '' );
    ref.a			= parameter( 'Linear Mic Array Angle', '\alpha', '°' );
    ref.T           = parameter( 'Jet Temperature', 'T', '°C' );
    ref.F           = parameter( 'Forcing Frequency','F','Hz');
    ref.S           = parameter( 'Forcing Strouhal Number','St_{DF}','');
    ref.x           = parameter( 'x/D of first microphone','x/D','');
    ref.r           = parameter( 'r/D of first microphone','r/D','');
    ref.PW          = parameter( 'Pulse width','PW','usec');
    ref.mV          = parameter( 'Calibration gains of microphones','mV','mV');
    ref.mVu         = parameter( 'Calibration gains of upstream microphones', 'mV','mV');
    ref.mVd         = parameter( 'Calibration gains of dowstream microphones', 'mV','mV');
    ref.mVf         = parameter( 'Calibration gains for farfield microphones', 'mV','mV');
end

function ref = NearFieldParams_v2()
    %Reference parameters   
    ref.M			= parameter( 'Mach Number', 'M', '' );
    ref.m           = parameter( 'Forcing Mode','m', '' );
    ref.a			= parameter( 'Linear Mic Array Angle', '\alpha', '°' );
    ref.T           = parameter( 'Jet Temperature', 'T', '°C' );
    ref.F           = parameter( 'Forcing Frequency','F','Hz');
    ref.S           = parameter( 'Forcing Strouhal Number','St_{DF}','');
    ref.x           = parameter( 'x/D of first microphone','x/D','');
    ref.r           = parameter( 'r/D of first microphone','r/D','');
    ref.PW          = parameter( 'Pulse width','PW','usec');
    ref.nca         = parameter( 'Calibration gains of microphones for first Nexus conditioner','mV','mV');
    ref.ncb         = parameter( 'Calibration gains of microphones for second Nexus conditioner','mV','mV');
    ref.ncc         = parameter( 'Calibration gains of microphones for third Nexus conditioner','mV','mV');
    ref.ncd         = parameter( 'Calibration gains of microphones for fourth Nexus conditioner','mV','mV');
    ref.nce         = parameter( 'Calibration gains of microphones for fifth Nexus conditioner','mV','mV');
end

function ref = MicCalibrationParams()
    %Reference parameters
    ref.dB          = parameter( 'Decibel Reference','dB','dB' );
    ref.Ch          = parameter( 'Mic Channel','Ch','' );
    ref.mVPa        = parameter( 'Mic Gain','mV/Pa','');
    ref.T           = parameter( 'Temperature','T','°C' );
    ref.CAL         = parameter( 'File Identifier','CAL','');
end

function ref = InletParams()
    %Reference parameters

    % Critical (these should always be included)
    ref.a			= parameter( 'Compression Angle', '\alpha', '°' );
    ref.M			= parameter( 'Mach Number', 'M', '' );
    ref.m           = parameter( 'Forcing Mode','m', '' );
    ref.p			= parameter( 'Stagnation Pressure', 'p_o', '% of 300 psig' );
    ref.pamb		= parameter( 'Ambient Pressure', 'p_{amb}', 'mbar' );
    ref.Tamb		= parameter( 'Ambient Temperature', 'T_{amb}', '°C' );

    % Non-critical
    ref.baseline	= parameter( 'Baseline', '', '' );
    ref.d			= parameter( 'Inviscid Impingement Point', 'd', 'mm' );
    ref.DC			= parameter( 'Duty Cycle', 'DC' ,'%' );
    ref.delta		= parameter( 'Boundary Layer Thickness', '\delta', 'mm' );
    ref.dt			= parameter( 'Laser Delay', 'dt', ' ns' );
    ref.Lint		= parameter( 'Interaction Length', 'L_{int}', 'mm' );
    ref.Lsep		= parameter( 'Separation Length', 'L_{sep}', 'mm' );
    ref.N			= parameter( 'Images in Ensemble', 'N', '' );
    ref.St			= parameter( 'Strouhal Number', 'St', '' );
    ref.T			= parameter( 'Stagnation Temperature', 'T_o', '°C' );
    ref.type		= parameter( 'PIV Type', '', '' );
    ref.Uinf		= parameter( 'Freestream Velocity', 'u_\infty', 'm/s' );

    % Special
    ref.U			= parameter( 'Streamwise Velocity', 'u', 'm/s' );
    ref.V			= parameter( 'Vertical Velocity', 'v', 'm/s' );
    ref.W			= parameter( 'Spanwise Velocity', 'w', 'm/s' );
    ref.X			= parameter( 'Streamwise Coordinate', 'x', 'mm' );
    ref.Y			= parameter( 'Vertical Coordinate', 'x', 'mm' );
    ref.Z			= parameter( 'Spanwise Coordinate', 'x', 'mm' );
end