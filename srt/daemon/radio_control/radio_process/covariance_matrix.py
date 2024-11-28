import numpy as np
from gnuradio import gr

class covariance_matrix_block(gr.sync_block):
    """
    A GNU Radio block that computes the covariance matrix elements
    for FFT streams across multiple channels.
    Each output corresponds to one element in the covariance matrix.
    """

    def __init__(self, num_channels=2, vec_length=256):
        """
        Initialize the block.
        Args:
            num_channels (int): Number of input channels.
            vec_length (int): Number of frequency bins per FFT.
        """
        gr.sync_block.__init__(
            self,
            name="Covariance Calculator Block",
            in_sig=[(np.complex64, vec_length)] * num_channels,
            out_sig=[(np.complex64, vec_length)] * (num_channels**2),
        )
        self.num_channels = num_channels
        self.vec_length = vec_length

    def work(self, input_items, output_items):
        """
        Perform the covariance matrix computation using full vectorization.
        """

        # Rearrange spectra into a form we can easily handle
        unwrapped_input = np.array(input_items).reshape(self.num_channels,-1)

        #combined_data = np.stack(input_items[:self.num_channels], axis=0)
        spectrum_matrices = np.transpose(unwrapped_input).reshape((-1, 1, self.num_channels))

        # Get covariance matrices (shape: vec_length, num_channels, num_channels)
        covariances = np.matmul(spectrum_matrices.swapaxes(1, 2), spectrum_matrices.conjugate())

        # The shape of covariances is now (vec_length, num_channels, num_channels)
        # We need to flatten it to the shape (num_channels * num_channels, vec_length)
        covariances_reshaped = (covariances.swapaxes(0, 2).reshape(self.num_channels**2, -1, self.vec_length))

        # Assign reshaped covariance matrix elements to the output items
        for idx in range(self.num_channels**2):
            output_items[idx][:] = covariances_reshaped[idx,:]

        # Return the number of items produced (number of bins)
        return len(output_items[0])
