function obj = select( obj, varargin )
%  SELECT - Select faces in COMPARTICLE object.
%
%  Usage for obj = comparticle :
%    obj = select( obj, 'PropertyName', PropertyValue )
%  PropertyName
%    'index'    :  index to selected elements
%    'carfun'   :  function f( x, y, z ) for selected elements
%    'polfun'   :  function f( phi, r, z ) for selected elements
%    'sphfun'   :  function f( phi, theta, r ) for selected elements
%  Output
%    obj        :  selected comparticle

if ~strcmp( varargin{ 1 }, 'index' )
  %  pass select input to all particle objects
  obj.p = cellfun( @( p ) select( p, varargin{ : } ), obj.p, 'UniformOutput', false );
else
  %  index to grouped particles
  ipt = cellfun( @( p, i ) repmat( i, 1, p.n ),  ...
                      obj.p, num2cell( 1 : numel( obj.p ) ), 'UniformOutput', false );
  ipt = horzcat( ipt{ : } );
  %  point index
  ind = cellfun( @( p ) 1 : p.n, obj.p, 'UniformOutput', false );  
  ind = horzcat( ind{ : } ); 
  %  get selected indices
  ind = ind( varargin{ 2 } );
  ipt = ipt( varargin{ 2 } );
  %  loop over all points
  for i = 1 : numel( obj.p )
    obj.p{ i } = select( obj.p{ i }, 'index', ind( ipt == i ) );
  end
end

%  index to non-empty particle objects
ind = cellfun( @( p ) ~isempty( p.pos ), obj.p );
%  keep only selected objects
[ obj.p, obj.inout ] = deal( obj.p( ind ), obj.inout( ind, : ) );
% %  reset closed arguments
obj.closed = cellfun( @( closed ) [], obj.closed, 'UniformOutput', false );
%  set mask object
obj.mask = 1 : numel( obj.p );

%  compound of particles
obj.pc = vertcat( obj.p{ obj.mask } );
