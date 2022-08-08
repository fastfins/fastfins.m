function f = eval_energy_rom(rom, Mbnd, fstate, state)

% energy norm
% f = 0;
f1 = (1/model.p_rate)*sum( ( fstate(:).^(model.p_rate/2) ).*model.WdetJ(:) );

% forcing
f2 = - model.b'*state;

f3 = 0.5 * state'*(Mbnd*state);

f = f1 + f2 + f3;

end