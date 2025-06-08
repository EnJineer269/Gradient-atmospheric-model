function T_interp = getAtmTransmittance(location,altitude,wavelengths)

persistent data originalNames;

% Load data once per session

if isempty(data)
    loaded = load('Atmospheric_Transmittance_Library.mat', 'atmTransmittanceData', 'originalNames');
    data = loaded.atmTransmittanceData;
    originalNames = loaded.originalNames;
end

% Make valid field names
locField = matlab.lang.makeValidName(location);
altField = ['alt_' strrep(num2str(altitude), '.', '_')];

% Check if location exists
if ~isfield(data, locField)
    error('Location "%s" not found in the library.', location);
end

% Check if altitude exists
if ~isfield(data.(locField), altField)
    error('Altitude %.1f km not found for location "%s".', altitude, location);
end

% Extract and interpolate
dataset = data.(locField).(altField);
T_interp = interp1(dataset.wavelength, dataset.transmittance, wavelengths, 'linear', 'extrap');

end


