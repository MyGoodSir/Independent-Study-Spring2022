"use strict";

class MeanCurvatureFlow {
	/**
	 * This class performs {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf mean curvature flow} on a surface mesh.
	 * @constructor module:Projects.MeanCurvatureFlow
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);
	}

	/**
	 * Builds the mean curvature flow operator.
	 * @private
	 * @method module:Projects.MeanCurvatureFlow#buildFlowOperator
	 * @param {module:LinearAlgebra.SparseMatrix} M The mass matrix of the input mesh.
	 * @param {number} h The timestep.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	buildFlowOperator(M, h) {
		return M.plus(this.geometry.laplaceMatrix(this.vertexIndex).timesReal(h));
	}

	/**
	 * Performs mean curvature flow on the input mesh with timestep h.
	 * @method module:Projects.MeanCurvatureFlow#integrate
	 * @param {number} h The timestep.
	 */
	integrate(h) {
		let vertices = this.geometry.mesh.vertices;
		let mass = this.geometry.massMatrix(this.vertexIndex);
		let op = this.buildFlowOperator(mass, h);
		let f=DenseMatrix.ones(vertices.length,3);
		for(let v of vertices){
			let p = this.geometry.positions[v], i = this.vertexIndex[v];
			f.set(p.x,i,0);
			f.set(p.y,i,1);
			f.set(p.z,i,2);
		}
		let A = op.chol();
		let b = mass.timesDense(f);
		let x = A.solvePositiveDefinite(b);
		for(let v of vertices){
			let i = this.vertexIndex[v];
			this.geometry.positions[v].x = x.get(i,0);
			this.geometry.positions[v].y = x.get(i,1);
			this.geometry.positions[v].z = x.get(i,2);
		}

		// center mesh positions around origin
		normalize(this.geometry.positions, vertices, false);
	}
}

