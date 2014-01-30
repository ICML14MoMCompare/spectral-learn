function [wa,A,n] = load_model(fname)
	% Get alphabet size and number of states
	fid = fopen(fname);
	n = 0;
	r = 0;
	[str,c] = fgets(fid);
	while (c > 0)
		% Parse string for initial and stopping probabilities
		[st,en] = regexp(str,'\([0-9]+\)');
		if (length(st) > 0)
			foo = sscanf(str(st(1):en(1)),'(%d)');
			n = max(n,foo);
		end;
		% Parse string for emission probabilities
		[st,en] = regexp(str,'\([0-9]+,[0-9]+\)');
		if (length(st) > 0)
			foo = sscanf(str(st(1):en(1)),'(%d,%d)');
			n = max(n,foo(1));
			r = max(r,foo(2));
		end;
		% Parse string for transition probabilities
		[st,en] = regexp(str,'\([0-9]+,[0-9]+,[0-9]+\)');
		if (length(st) > 0)
			foo = sscanf(str(st(1):en(1)),'(%d,%d,%d)');
			n = max(n,foo(1));
			r = max(r,foo(2));
			n = max(n,foo(3));
		end;
		% Read next string
		[str,c] = fgets(fid);
	end;
	fclose(fid);
	% XXX Beware! Symbols and states are indexed from 0
	n = n + 1;
	r = r + 1;
	% Build alphabet
	A = [];
	for i = 1:r
		A = [A char('a' + i - 1)];
	end;
	% Read file again, build matrices and tensors
	a1 = zeros(1,n);
	ainf = zeros(n,1);
	O = zeros(n,r);
	T = zeros(n,r,n);
	fid = fopen(fname);
	[str,c] = fgets(fid);
	initial = true;
	while (c > 0)
		% Switch from initial to stopping probabilities
		[st,en] = regexp(str,'F:');
		if (length(st) > 0)
			initial = false;
		end;
		% Parse string for initial and stopping probabilities
		[st,en] = regexp(str,'\([0-9]+\)');
		if (length(st) > 0)
			foo = sscanf(str(st(1):en(1)),'(%d)');
			bar = sscanf(str(en(1)+1:end),'%f');
			if (initial)
				a1(foo+1) = bar;
			else
				ainf(foo+1) = bar;
			end;
		end;
		% Parse string for emission probabilities
		[st,en] = regexp(str,'\([0-9]+,[0-9]+\)');
		if (length(st) > 0)
			foo = sscanf(str(st(1):en(1)),'(%d,%d)');
			bar = sscanf(str(en(1)+1:end),'%f');
			O(foo(1)+1,foo(2)+1) = bar;
		end;
		% Parse string for transition probabilities
		[st,en] = regexp(str,'\([0-9]+,[0-9]+,[0-9]+\)');
		if (length(st) > 0)
			foo = sscanf(str(st(1):en(1)),'(%d,%d,%d)');
			bar = sscanf(str(en(1)+1:end),'%f');
			T(foo(1)+1,foo(2)+1,foo(3)+1) = bar;
		end;
		% Read next string
		[str,c] = fgets(fid);
	end;
	fclose(fid);
	% Normalize O using stopping probabilities
	N = zeros(size(O));
	for i = 1:size(O,2)
		N(:,i) = 1 - ainf;
	end;
	O = O .* N;
	% Convert O and T to WFA parametrization
	wa = struct();
	wa.a1 = a1;
	wa.ainf = ainf;
	for k = 1:length(A)
		a = A(k);
		M = zeros(n,n);
		for i = 1:n
			for j = 1:n
				M(i,j) = O(i,k) * T(i,k,j);
			end;
		end;
		wa.(a) = M;
	end;


