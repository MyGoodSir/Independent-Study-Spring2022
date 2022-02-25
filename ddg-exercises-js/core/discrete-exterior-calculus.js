"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {
		let trp = new Triplet(geometry.mesh.vertices.length, geometry.mesh.vertices.length);
		geometry.mesh.vertices.forEach(v => {
			trp.addEntry(
				geometry.barycentricDualArea(v),
				vertexIndex[v],
				vertexIndex[v]
			);
		});
		return SparseMatrix.fromTriplet(trp);
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
		let trp = new Triplet(geometry.mesh.edges.length, geometry.mesh.edges.length);
		geometry.mesh.edges.forEach(e => {
			trp.addEntry(
				(geometry.cotan(e.halfedge)+geometry.cotan(e.halfedge.twin))/2.0,
				edgeIndex[e],
				edgeIndex[e]
			);
		});
		return SparseMatrix.fromTriplet(trp);
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
		let trp = new Triplet(geometry.mesh.faces.length, geometry.mesh.faces.length);
		geometry.mesh.faces.forEach(f => {
			trp.addEntry(
				(a => {
					let hsum=a.reduce((x,b)=>x+b)/2;
					return 1/Math.sqrt(a.reduce((acc, l)=>{return acc*(hsum-l)},hsum))
				})([...f.adjacentEdges()].map(e=>geometry.length(e))),
				faceIndex[f],
				faceIndex[f]
			);
		});
		return SparseMatrix.fromTriplet(trp);
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		let trp = new Triplet(geometry.mesh.edges.length, geometry.mesh.vertices.length)
		geometry.mesh.vertices.forEach(v=>{[...v.adjacentEdges()].forEach(e=>trp.addEntry(-2*(vertexIndex[e.halfedge.vertex]==vertexIndex[v])+1,edgeIndex[e],vertexIndex[v]))})
		return SparseMatrix.fromTriplet(trp);
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
		let trp = new Triplet(geometry.mesh.faces.length, geometry.mesh.edges.length)
		geometry.mesh.faces.forEach(f=>{[...f.adjacentEdges()].forEach(e=>trp.addEntry(2*(faceIndex[e.halfedge.face]==faceIndex[f])-1,faceIndex[f],edgeIndex[e]))})
		return SparseMatrix.fromTriplet(trp);
	}
}
