import numpy as np
from gnuradio import gr

class covariance_matrix_block(gr.sync_block):
    """
    A GNU Radio block that computes the covariance matrix elements
    for FFT streams across multiple channels.
    Each output corresponds to one element in the covariance matrix.
    """

    def __init__(self, num_channels=2, num_bins=256):
        """
        Initialize the block.
        Args:
            num_channels (int): Number of input channels.
            num_bins (int): Number of frequency bins per FFT.
        """
        gr.sync_block.__init__(
            self,
            name="Covariance Calculator Block",
            in_sig=[(np.complex64, num_bins)] * num_channels,
            out_sig=[(np.complex64, num_bins)] * (num_channels * num_channels),
        )
        self.num_channels = num_channels
        self.num_bins = num_bins

    def work(self, input_items, output_items):
        """
        Perform the covariance matrix computation using full vectorization.
        """

        # Rearrange spectra into a form we can easily handle
        spectrum_matrices = np.transpose(input_items).reshape((-1, 1, self.num_channels))

        # Get covariance matrices (shape: num_bins, num_channels, num_channels)
        covariances = np.matmul(spectrum_matrices.swapaxes(1, 2), spectrum_matrices.conjugate())

        # The shape of covariances is now (num_bins, num_channels, num_channels)
        # We need to flatten it to the shape (num_channels * num_channels, num_bins)
        covariances_reshaped = covariances.swapaxes(0, 2).reshape(self.num_channels**2, -1)

        # Assign reshaped covariance matrix elements to the output items
        for idx in range(self.num_channels**2):
            output_items[idx][:] = covariances_reshaped[idx, :]

        # Return the number of items produced (number of bins)
        return len(output_items[0])
