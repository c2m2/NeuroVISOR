﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using TMPro;

namespace C2M2.NeuronalDynamics.Visualization
{
    using Interaction;
    using VRN;
    using Utils.DebugUtils;
    using Utils;
    using UGX;
    using Grid = UGX.Grid;
    using DiameterAttachment = UGX.IAttachment<UGX.DiameterData>;

    /// <summary>
    /// Produces a preview of a 1D cell using LinesRenderer
    /// </summary>
    public class NeuronCellPreview : MonoBehaviour
    {
        public string vrnFileName = "null";
        public Color32 color;
        public NDSimulationLoader loader = null;
        public TextMeshProUGUI fileNameDisplay;
        public TextMeshProUGUI sizeLabel;
        public TextMeshProUGUI[] textColTargets = new TextMeshProUGUI[0];
        public string LengthScale { get { return loader.lengthScale; } }
        public int refinement = 0;
        public int[] refinements { get; private set; }
        private VrnReader vrnReader = null;
        private LinesRenderer lines = null;

        public void PreviewCell()
        {
            PreviewCell(vrnFileName, color);
        }
        public void PreviewCell(string vrnFileName, Color32 color)
        {
            if(vrnFileName == "null")
            {
                Debug.LogError("Null cell given to NeuronCellPreview");
                if (fileNameDisplay != null) fileNameDisplay.text = "null";
                return;
            }

            char sl = Path.DirectorySeparatorChar;
            if (!vrnFileName.EndsWith(".vrn")) vrnFileName = vrnFileName + ".vrn";
            vrnReader = new VrnReader(Application.streamingAssetsPath + sl + "NeuronalDynamics" + sl + "Geometries" + sl + vrnFileName);

            refinements = vrnReader.ListRefinements();

            string meshName1D = vrnReader.Retrieve1DMeshName(refinement);

            /// Create empty grid with name of grid in archive
            Grid grid = new Grid(new Mesh(), meshName1D);
            grid.Attach(new DiameterAttachment());

            // Read the cell
            vrnReader.ReadUGX(meshName1D, ref grid);

            // Scale the parent object by 1 / max scale to make the cell fit within size (1,1,1)
            float scale = 1 / Math.Max(grid.Mesh.bounds.size);
            transform.localScale = new Vector3(scale, scale, scale);

            // Adjust center so cell mesh is centered at (0,0,0)
            transform.localPosition = -scale * grid.Mesh.bounds.center;

            if(lines != null)
            {
                Destroy(lines);
                lines = null;
            }

            // Render cells
            lines = gameObject.AddComponent<LinesRenderer>();

            // (line width = scale)
            lines.Draw(grid, color, scale);

            // Displays file name without file extension
            if (fileNameDisplay != null)
                fileNameDisplay.text = vrnFileName.Remove(vrnFileName.LastIndexOf('.'));

            // Draw scale labels
            Vector3 cellSize = grid.Mesh.bounds.size;
            if (sizeLabel != null)
                sizeLabel.text = 
                    "Size: ("
                    + cellSize.x.ToString() + ", "
                    + cellSize.y.ToString() + ", " 
                    + cellSize.z.ToString() + " " + LengthScale + ")";

            foreach(var t in textColTargets)
            {
                t.color = color;
            }
        }
        public void LoadThisCell(RaycastHit hit)
        {
            loader.vrnFileName = vrnFileName;
            loader.refinementLevel = refinement;
            loader.Load(hit);
        }
        public bool RemoveRefinement(int refinement)
        {
            bool removeSuccessful = false;

            List<int> refOptions = refinements.ToList();
            if (refOptions.Contains(refinement))
            {
                refOptions.Remove(refinement);

                refinements = refOptions.ToArray();
                removeSuccessful = true;
            }

            if(this.refinement == refinement)
            {
                if(this.refinement-1 > 0)
                {
                    this.refinement--;
                    PreviewCell();
                } 
            }

            return removeSuccessful;
        }
    }
}