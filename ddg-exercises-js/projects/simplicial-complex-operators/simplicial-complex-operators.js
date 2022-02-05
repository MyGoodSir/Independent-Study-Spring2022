"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                mesh.edges.forEach((e, i) => e.index = i);
                mesh.vertices.forEach((v, i) => v.index = i);
                mesh.faces.forEach((f, i) => f.index = i);
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                let triplet = new Triplet(mesh.edges.length, mesh.vertices.length);
                mesh.edges.forEach(e => {
                        triplet.addEntry(1, e.index, e.halfedge.vertex.index); 
                        triplet.addEntry(1, e.index, e.halfedge.twin.vertex.index);
                });
                return SparseMatrix.fromTriplet(triplet);
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                let triplet = new Triplet(mesh.faces.length, mesh.edges.length);
                mesh.faces.forEach(f => [...f.adjacentEdges()].forEach(e => triplet.addEntry(1, f.index, e.index)));
                return SparseMatrix.fromTriplet(triplet);
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                return ((m) => {
                        subset.vertices.forEach(v=>m.set(1,v,0));return m;
                })(DenseMatrix.zeros(this.mesh.vertices.length, 1))
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                return ((m) => {
                        subset.edges.forEach(e=>m.set(1,e,0));return m;
                })(DenseMatrix.zeros(this.mesh.edges.length, 1))
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                return ((m) => {
                        subset.faces.forEach(f=>m.set(1,f,0));return m;
                })(DenseMatrix.zeros(this.mesh.faces.length, 1))
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                return ((s) => {
                        let e = this.buildEdgeVector(subset).plus(this.A0.timesDense(this.buildVertexVector(subset)));
                        [...e.values()].map((v, i) => {if(v!=0){s.addEdge(i)}});
                        [...this.buildFaceVector(subset).plus(this.A1.timesDense(e)).values()].map((v,i) => {if(v!=0){s.addFace(i)}});
                        return s;
                })(MeshSubset.deepCopy(subset));
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                return ((s) => {
                        let e = this.buildEdgeVector(subset).plus(this.A1.transpose().timesDense(this.buildFaceVector(subset)));
                        [...e.values()].map((v, i) => {if(v!=0){s.addEdge(i)}});
                        [...this.buildVertexVector(subset).plus(this.A0.transpose().timesDense(e)).values()].map((v,i) => {if(v!=0){s.addVertex(i)}});
                        return s;
                })(MeshSubset.deepCopy(subset));
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                let out = this.closure(this.star(subset));
                out.deleteSubset(this.star(this.closure(subset)));
                return out; // placeholder
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                return subset.equals(this.closure(subset));
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                if(!this.isComplex(subset)){ return -1;}

                let e = this.buildEdgeVector(subset), v = this.buildVertexVector(subset), degree = 0;

                if (subset.faces.size > 0 ){
                        let fe = this.A1.transpose().timesDense(this.buildFaceVector(subset));
                        if([...Array(e.nRows()).keys()].some(i => (fe.get(i,0)==0) != (e.get(i,0)==0))){ return -1;}
                        degree = 2;
                }

                let ev = this.A0.transpose().timesDense(e);

                if (subset.edges.size > 0){
                        if([...Array(v.nRows()).keys()].some(i => (ev.get(i,0)==0) != (v.get(i,0)==0))){ return -1;}
                        if (degree!= 2){ degree = 1;}
                }

                return degree;
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                let fn = subset.faces.size > 0 ? ((b)=>{
                        [...this.A1.transpose().timesDense(this.buildFaceVector(subset)).values()].forEach((v,i)=>{if(v==1){b.addEdge(i)}}); return b;
                }) : (subset.edges.size > 0 ? (b)=>{ 
                        [...this.A0.transpose().timesDense(this.buildEdgeVector(subset)).values()].forEach((v,i)=>{if(v==1){b.addVertex(i)}}); return b;
                } : (b)=>b)
                return this.closure(fn(new MeshSubset()));
        }
}
