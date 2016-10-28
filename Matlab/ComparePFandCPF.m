% Script to compare CPF and repetitive load flows

mpopt = mpoption('out.all', 0, 'verbose', 2);
mpopt = mpoption(mpopt, 'cpf.stop_at', 'FULL', 'cpf.step', 0.2);
mpopt = mpoption(mpopt,'cpf.plot.level',2,'cpf.parameterization',3);
runcpf('case9','case9target',mpopt);