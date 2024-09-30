# MCNET - Gene Regulatory Network Mapping

MCNET is a state-of-the-art deep learning model designed to infer gene regulatory networks (GRNs) using multi-omics integration from single-cell RNA sequencing (scRNA-seq) data. By incorporating attention mechanisms and graph convolutional networks, MCNET efficiently captures the intricate regulatory relationships among genes. It also utilizes copy number variations (CNVs) and DNA methylation data to enhance the accuracy of GRN inference.

## Key Features

- **Multi-Omics Integration**: MCNET integrates data from multiple sources, including scRNA-seq, CNV, and DNA methylation to improve gene regulatory network inference.
- **Deep Learning Architecture**: The model leverages graph convolutional networks and attention mechanisms to better capture complex gene interactions.
- **GRN Inference**: Accurately predicts gene regulations on cell-type marker genes using multi-layered adjacency matrices.
- **Visualization & Clustering**: Includes robust mechanisms for visualizing scRNA-seq data and clustering based on cell-type identification.
- **Simulation Studies**: Includes simulation data to validate and benchmark the performance of the model against state-of-the-art techniques.

## Technology Stack

- **Deep Learning**: Implemented using TensorFlow/PyTorch frameworks for neural network computations.
- **Graph Convolutional Networks (GCN)**: For modeling gene interactions and constructing the gene regulatory networks.
- **Attention Mechanisms**: Incorporated to focus on critical interactions within the multi-omics data.
- **Multi-Omics Data Handling**: Supports integration of multiple omics datasets including scRNA-seq, CNVs, and DNA methylation.
- **Python**: Core language used for development.
- **Matplotlib/Seaborn**: For data visualization and generating figures.

## Dataset

MCNET was tested and validated on multiple publicly available single-cell RNA sequencing datasets, along with CNV and DNA methylation data:

- **scRNA-seq datasets**: Mouse cortex cell-type marker genes.
- **Epigenetic data**: Validates the gene regulatory predictions.

## Setup and Installation

### Prerequisites

- Python 3.x
- TensorFlow or PyTorch (based on the chosen deep learning framework)
- Required Python libraries: `numpy`, `pandas`, `scikit-learn`, `matplotlib`, `seaborn`

### Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/ansschh/MCNET-Gene_Regulatory_Network_Mapping.git
    cd MCNET-Gene_Regulatory_Network_Mapping
    ```

2. Install the required Python packages:
    ```bash
    pip install -r requirements.txt
    ```

3. Download or prepare your dataset (scRNA-seq, CNV, and DNA methylation data).

### Usage

To train the MCNET model and infer gene regulatory networks:

1. **Prepare the dataset**: Ensure your data is in the correct format (e.g., gene expression matrix, CNV data, DNA methylation data).

2. **Run the training script**:
    ```bash
    python train_mcnet.py --data_path /path/to/dataset --epochs 50
    ```

3. **Visualize the results**: Once the model is trained, use the provided visualization tools to interpret the inferred gene regulatory networks.
    ```bash
    python visualize_grn.py --output_path /path/to/save/figures
    ```

4. **Benchmark with Simulation Data**: Run simulation experiments to evaluate the model’s performance.
    ```bash
    python simulate.py --simulation_type roc_curve
    ```
To run the program with the introduced data in the article, you only need to do as follows:

1- Run R 

2- library(lqa) 

3- trace(fused.lasso,edit=T)

4- substitute the existing function with the function in fused.lasso.modified.r (following function):
						
	function (lambda = NULL, ...) 					
	{					
	  argList = list(...)					
	  w <- argList$...					
	  lambda.check(lambda)					
	  if (length(lambda) != 2) 					
	    stop("The fused.lasso penalty must consist on two parameters! \n")					
	  names(lambda) <- c("lambda1", "lambda2")					
	  first.derivative <- function(beta = NULL, ...) {					
	    if (is.null(beta)) 					
	      stop("'beta' must be the current coefficient vector \n")					
	    p <- length(beta)					
	    if (p < 2) 					
	      stop("There must be at least two regressors! \n")					
	    vec1 <- c(rep(lambda[1], p), rep(lambda[2], (3/4) * p))					
	    vec2 <- abs(drop(t(a.coefs(beta, ... = w)) %*% beta))					
	    #print(any(vec2 < 0))					
	    return(vec1 * vec2)					
	  }					
	  a.coefs <- function(beta = NULL, ...) {					
	    argList = list(...)					
	    w <- argList$...					
	    if (is.null(beta)) 					
	      stop("'beta' must be the current coefficient vector \n")					
	    p <- length(beta)					
	    if (p < 2) 					
	      stop("There must be at least two regressors! \n")					
	    if (p > 2) {					
	      h1 <- cbind(-diag((3/4) * p), matrix(0, (3/4) * p, (1/4) * p))					
	      h2 <- cbind(matrix(0, (3/4) * p, (1/4) * p), diag((3/4) * p))					
	      #print("i did it :)")					
	      mat1 <- h1 + h2					
	      mat2 <- diag(w)					
	      a.coefs.mat <- cbind(mat2, t(mat1))					
	    }					
	    else a.coefs.mat <- cbind(diag(2), c(-1, 1))					
	    return(a.coefs.mat)					
	  }					
	  structure(list(penalty = "fused.lasso", lambda = lambda, 					
		         first.derivative = first.derivative, a.coefs = a.coefs), 				
		    class = "penalty")				
	}					

5- source("sources.r") #make sure that source.r is in the correct directory



## Results

- **GRN Inference Accuracy**: MCNET outperforms existing methods with an AUROC score of 0.938 compared to GRN (0.895) and DCGRN (0.843).
- **True Positive Rate (TPR)**: MCNET shows superior performance in identifying true positive interactions from multi-omics data.
- **Sensitivity**: MCNET achieves better sensitivity (0.300 ± 0.034) in simulations compared to baseline methods.
- **False Discovery Rate (FDR)**: Maintains low FDR (<0.02) across experiments, indicating reliable gene regulatory network predictions.

## Citation

If you use MCNET in your research, please cite the following paper:

_Tiwari, A., & Trankatwar, S. (2023). MCNET: Multi-Omics Integration for Gene Regulatory Network Inference from scRNA-seq. Birla Institute of Technology & Science Pilani, Hyderabad Campus._


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

We welcome contributions to enhance MCNET! Please follow the steps below to contribute:

1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/new-feature`).
3. Commit your changes (`git commit -m 'Add some feature'`).
4. Push to the branch (`git push origin feature/new-feature`).
5. Open a pull request.

## Contact

For any questions or issues, feel free to open an issue or contact the project maintainers:

- Ansh Tiwari (anshtiwari9899@gmail.com)
- Sachin Trankatwar (strankatwar@gmail.com)
